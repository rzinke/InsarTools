#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Extract profiles from a MintPy displacement data cube.

FUTURE IMPROVEMENTS
    * Plot profiles

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
from IOsupport import load_profile_endpoints, confirm_outdir, confirm_outname_ext, append_fname, load_mintpy_timeseries, save_profile_data
from GeoFormatting import get_mintpy_transform, transform_to_extent, lola_to_xy, get_mintpy_reference_point
from Masking import create_mask
from Profiling import profile_geometry, format_profile_grid, extract_profile_values
from Fitting import dates_to_datetimes, time_since_reference, fit_atan
from Viewing import plot_raster, plot_profile


### PARSER ---
Description = '''Extract profiles from a MintPy displacement data cube.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='tsName', type=str,
        help='MintPy timeseries.h5 file.')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')


    QueryArgs = parser.add_argument_group('QUERY ARGUMENTS')
    QueryArgs.add_argument('-q','--queryfile', dest='queryFile', type=str, default=None,
        help='Profile endpoints in geographic coordinates, one profile per line (lonStart,latStart lonEnd,latEnd).')
    QueryArgs.add_argument('-w','--profile-width', dest='profWidth', type=float, default=None,
        help='Profile width in map units.')
    QueryArgs.add_argument('--fit-type', dest='fitType', type=str, default=None,
        help='Fit function ([None], linear, seasonal).')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true', 
        help='Plot data.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### TS PROFILE CLASS ---
class MintPyTSprofile:
    def __init__(self, tsName, outName, maskArgs=None, fitType=None, verbose=False):
        '''
        Extract 1D profiles from a MintPy timeseries.h5 data cube.
        '''
        # Parameters
        self.verbose = verbose

        # Format fit type
        self.__format_fit_type__(fitType)

        # Load MintPy data set
        self.__load_data__(tsName)

        # Create mask
        self.mask = create_mask(self.disps[-1,:,:], maskArgs, verbose=self.verbose)

        # Retrieve spatial extent
        self.__parse_spatial_info__(tsName)

        # Build grid
        self.__build_grid__()


    def __load_data__(self, tsName):
        '''
        Load the MintPy (HDF5) data set, parse the dates and displacements.
        '''
        # Open and close HDF5 data set
        self.dates, self.disps = load_mintpy_timeseries(tsName, verbose=self.verbose)

        # Convert dates to datetimes
        self.datetimes = dates_to_datetimes(self.dates, verbose=self.verbose)

        # Calculate time since beginning of series
        self.times = np.array(time_since_reference(self.datetimes, verbose=self.verbose))


    def __format_fit_type__(self, fitType):
        '''
        Check and format the type of fit to apply, if applicable.
        '''
        if fitType is not None:
            # Ensure lower case
            fitType = fitType.lower()

            # Check that fit type is acceptable
            assert fitType in ['atan', 'arctan'], \
                'Fit type {:s} is not valid, use \'linear\' or \'seasonal\''.format(fitType)

            # Confirm formatting
            if fitType in ['atan', 'arctan']: fitType = 'atan'

        self.fitType = fitType


    def __parse_spatial_info__(self, tsName):
        '''
        Retrieve the geographic information
        '''
        # Data set sizes
        self.nEpochs, self.M, self.N = self.disps.shape

        # Geographic transform
        self.tnsf = get_mintpy_transform(tsName, verbose=inps.verbose)

        # Geographic extent
        self.extent = transform_to_extent(self.tnsf, self.M, self.N, verbose=self.verbose)

        # Reference point
        self.refLon, self.refLat = get_mintpy_reference_point(tsName, verbose=self.verbose)


    def __build_grid__(self):
        '''
        Build an X, Y grid based on the data set size.
        '''
        x = np.arange(self.N)
        y = np.arange(self.M)

        self.X, self.Y = np.meshgrid(x, y)


    ## Query profiles
    def query_profiles(self, queryFile, profWidth):
        '''
        Query the profile between the points given in qLoLa and with the
         specified width.
        '''
        if self.verbose == True: print('Querying location:')

        # Parse profile points
        self.__parse_prof_points__(queryFile)

        # Format grid
        self.__format_grid__()

        # Format profile width
        self.__format_width__(profWidth)

        # Compute full profile geometry
        self.profGeom = profile_geometry(self.verbose)
        self.profGeom.from_endpoints((self.startLon, self.startLat), (self.endLon, self.endLat), self.profWidth)

        # Loop through each epoch to extract profile points
        self.profDists = []
        self.profPts = []
        pxSize = self.tnsf[1]

        for i in range(self.nEpochs):
            if self.verbose == True: print('Extracting profile for {:s}'.format(self.dates[i]))

            pxDist, profPts = extract_profile_values(img=self.disps[i,:,:],
                X=self.X, Y=self.Y,
                pxLen=self.pxLen, pxWidth=self.pxWidth,
                mask=self.mask,
                verbose=self.verbose)

            # Scale profile distance from pixels to map units
            profDist = pxDist*pxSize

            # Append to lists
            self.profDists.append(profDist)
            self.profPts.append(profPts)


    def __parse_prof_points__(self, queryFile):
        '''
        Confirm that the query points are in the correct format, i.e.,
            ((startLon,startLat) (endLon,endLat))
        and are within the extent of the map.
        '''
        # Load profile inputs
        startLons, startLats, endLons, endLats = load_profile_endpoints(queryFile, verbose=self.verbose)

        if len(startLons) > 1: print('WARNING: More than one profile detected.')

        # Parse profile inputs
        self.startLon = startLons[0]
        self.startLat = startLats[0]

        self.endLon = endLons[0]
        self.endLat = endLats[0]

        # Report if requested
        if self.verbose == True:
            print('Profile start: {:f}, {:f}'.format(self.startLon, self.startLat))
            print('Profile end: {:f}, {:f}'.format(self.endLon, self.endLat))

        # Convert to image coordinates
        self.pxStart, self.pyStart = lola_to_xy(self.tnsf, self.startLon, self.startLat,
            verbose=self.verbose)
        self.pxEnd, self.pyEnd = lola_to_xy(self.tnsf, self.endLon, self.endLat,
            verbose=self.verbose)

        # Check that queried coordinates are within bounds
        assert min([self.pxStart, self.pyStart, self.pxEnd, self.pyEnd]) >= 0, \
            'All coordinates must be within map bounds'

        assert max([self.pxStart, self.pxEnd]) <= self.N, 'X-coordinates must be within map bounds'
        assert max([self.pyStart, self.pyEnd]) <= self.M, 'Y-coordinates must be within map bounds'


    def __format_grid__(self):
        '''
        Format the pixel grid such that X, Y are centered and rotated such that
         X increases with distance along the profile length, and Y increases or
         decreases with profile width.
        Do this here so that it only needs to be done once per data set.
        '''
        self.X, self.Y, self.pxLen = format_profile_grid(self.X, self.Y,
            self.pxStart, self.pyStart, self.pxEnd, self.pyEnd,
            verbose=self.verbose)


    def __format_width__(self, profWidth):
        '''
        Format the profile with based on the pixel size if no set width is
         specified.
        '''
        # Pixel size
        pxSize = self.tnsf[1]

        # If width is not specified, use intrinsic pixel size
        if profWidth is None:
            self.profWidth = 2*pxSize  # use x-step
        else:
            self.profWidth = profWidth

        # Convert to pixels
        self.pxWidth = int(self.profWidth/pxSize)


    ## Plotting
    def plot(self):
        '''
        Plot profile locations and profile data.
        '''
        # Plot basemap
        mapFig, axMap = plot_raster(self.disps[-1,:,:],
            mask=self.mask, extent=self.extent,
            minPct=1, maxPct=99, cbarOrient='vertical')

        # Plot reference point
        axMap.plot(self.refLon, self.refLat, 'ks')

        # Plot profile geometry
        plot_profile(self.profGeom, fig=mapFig, ax=axMap)
        axMap.plot(self.startLon, self.startLat, 'ks')
        axMap.plot([self.startLon, self.endLon], [self.startLat, self.endLat], 'k')


    ## Saving
    def save(self, outName):
        '''
        Save each profile to a text file.
        '''
        if self.verbose == True: print('Saving to files using prefix {:s}'.format(outName))

        # Confirm output directory exists
        confirm_outdir(outName)

        # Confirm extention
        outName = confirm_outname_ext(outName, ['txt'], self.verbose)

        # Loop through profiles
        for i in range(self.nEpochs):
            # Append profile name
            fname = append_fname(outName, '_{:s}'.format(self.dates[i]))

            # Report if requested
            if self.verbose == True: print('... {:s}'.format(os.path.basename(fname)))

            # Save to file
            save_profile_data(fname, self.profGeom.profStart, self.profGeom.profEnd,
                self.profDists[i], self.profPts[i])



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Extract profiles
    TS = MintPyTSprofile(tsName=inps.tsName, outName=inps.outName, maskArgs=inps.maskArgs, fitType=inps.fitType,
        verbose=inps.verbose)

    # Query profiles
    TS.query_profiles(inps.queryFile, inps.profWidth)

    # Plot data
    if inps.plot == True:
        TS.plot()

    # Save to files
    if inps.outName is not None:
        TS.save(inps.outName)


    plt.show()