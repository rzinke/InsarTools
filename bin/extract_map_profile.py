#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Retrieve a profile from a map.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import load_profile_endpoints, load_gdal_dataset, confirm_outdir, confirm_outname_ext, confirm_overwrite, save_profile_data
from GeoFormatting import transform_to_extent, lola_to_xy
from Masking import create_mask
from Profiling import profile_geometry, extract_profile
from Viewing import plot_raster, plot_profile


### PARSER ---
Description = '''Examine a profile along a map.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='imgName', type=str,
        help='Image file name.')
    InputArgs.add_argument('-b','--band', dest='band', type=int, default=1,
        help='Image band number ([1]).')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')

    ProfileArgs = parser.add_argument_group('PROFILE ARGUMENTS')
    ProfileArgs.add_argument('-q','--queryfile', dest='queryFile', type=str, default=None,
        help='Profile endpoints in geographic coordinates, one profile per line (lonStart,latStart lonEnd,latEnd).')
    ProfileArgs.add_argument('-w','--profile-width', dest='profWidth', type=float, default=None,
        help='Profile width in map units.')

    DisplayArgs = parser.add_argument_group('DISPLAY ARGUMENTS')
    DisplayArgs.add_argument('-c','--cmap', dest='cmap', type=str, default='viridis',
        help='Colormap')
    DisplayArgs.add_argument('-co','--cbar-orient', dest='cbarOrient', type=str,
        default='horizontal', help='Colorbar orientation ([horizontal], vertical')
    DisplayArgs.add_argument('-minPct','--min-percent', dest='minPct', type=float, default=None,
        help='Minimum percent clip value ([None]).')
    DisplayArgs.add_argument('-maxPct','--max-percent', dest='maxPct', type=float, default=None,
        help='Maximum percent clip value ([None]).')

    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true', 
        help='Plot map and outputs.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default=None, 
        help='Output name.')
    OutputArgs.add_argument('--overwrite', dest='overwrite', action='store_true',
        help='Force overwrite output file.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### MAP PROFILE CLASS ---
class mapProfile:
    def __init__(self, imgName, band=1, maskArgs=None, verbose=False):
        '''
        Load an image and collect a profile across it.
        '''
        # Parameters
        self.verbose = verbose

        # Load data set
        self.__load_data__(imgName, band)

        # Create mask
        self.mask = create_mask(self.img, maskArgs, verbose=self.verbose)


    ## Loading
    def __load_data__(self, imgName, band):
        '''
        Load image data set using GDAL and retrieve geospatial information.
        '''
        if self.verbose == True: print('Loading: {:s}'.format(imgName))

        # Load data set using GDAL
        DS = load_gdal_dataset(imgName)

        # Extract image values
        self.img = DS.GetRasterBand(band).ReadAsArray()

        # Replace NaNs with zeros
        self.img[np.isnan(self.img) == 1] = 0

        # Extract geographic information
        self.M, self.N = DS.RasterYSize, DS.RasterXSize
        self.tnsf = DS.GetGeoTransform()
        self.extent = transform_to_extent(self.tnsf, self.M, self.N, verbose=self.verbose)        


    ## Profiling
    def query_profile(self, queryFile, profWidth):
        '''
        Query the profile between the points given in qLoLa and with the
         specified width.
        '''
        if self.verbose == True: print('Querying location:')

        # Parse profile points
        self.__parse_prof_points__(queryFile)

        # Format profile width
        self.__format_width__(profWidth)

        # Compute full profile geometry
        self.profGeom = profile_geometry(self.verbose)
        self.profGeom.from_endpoints((self.startLon, self.startLat), (self.endLon, self.endLat), self.profWidth)

        # Extract profile
        self.pxDist, self.profPts = extract_profile(img=self.img,
            pxStart=self.pxStart, pyStart=self.pyStart,
            pxEnd=self.pxEnd, pyEnd=self.pyEnd,
            pxWidth=self.pxWidth,
            mask=self.mask,
            verbose=self.verbose)

        # Scale profile distance from pixels to map units
        pxSize = self.tnsf[1]
        self.profDist = self.pxDist*pxSize

        # Number of profile points
        self.nPts = len(self.profPts)


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
        pxStart, pyStart = lola_to_xy(self.tnsf, np.array(self.startLon), np.array(self.startLat),
            verbose=self.verbose)
        pxEnd, pyEnd = lola_to_xy(self.tnsf, np.array(self.endLon), np.array(self.endLat),
            verbose=self.verbose)

        self.pxStart = pxStart[0]; self.pyStart = pyStart[0]
        self.pxEnd = pxEnd[0]; self.pyEnd = pyEnd[0]

        # Check that queried coordinates are within bounds
        assert min([self.pxStart, self.pyStart, self.pxEnd, self.pyEnd]) >= 0, \
            'All coordinates must be within map bounds'

        assert max([self.pxStart, self.pxEnd]) <= self.N, 'X-coordinates must be within map bounds'
        assert max([self.pyStart, self.pyEnd]) <= self.M, 'Y-coordinates must be within map bounds'


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
    def plot(self, cmap, cbarOrient, minPct, maxPct):
        '''
        Plot the map and profile.
        '''
        # Create map plot
        mapFig, axMap = plt.subplots(figsize=(8,8))

        # Plot map
        mapFig, axMap = plot_raster(self.img, mask=self.mask,
            extent=self.extent,
            cmap=cmap, cbarOrient=cbarOrient,
            minPct=minPct, maxPct=maxPct,
            fig=mapFig, ax=axMap)

        # Plot profile
        plot_profile(self.profGeom, fig=mapFig, ax=axMap)
        axMap.plot(self.startLon, self.startLat, 'ks')
        axMap.plot([self.startLon, self.endLon], [self.startLat, self.endLat], 'k')

        # Create profile plot
        profFig, axProf = plt.subplots(figsize=(9,5))
        axProf.plot(self.profDist, self.profPts, 'k.')
        axProf.set_xlabel('distance along profile')


    ## Save to file
    def save(self, outName, overwrite=False):
        '''
        Save profile points to file.
        '''
        if self.verbose == True: print('Saving to file: {:s}'.format(outName))

        # Confirm output directory exists
        confirm_outdir(outName)

        # Confirm extention
        outName = confirm_outname_ext(outName, ['txt'], self.verbose)

        # Check overwrite status
        if overwrite == False:
            overwrite = confirm_overwrite(outName)

        if overwrite == True:
            save_profile_data(outName, (self.startLon, self.startLat), (self.endLon, self.endLat),
                self.profDist, self.profPts, self.verbose)
        else:
            print('Saving aborted.')



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Extract profile
    prof = mapProfile(imgName=inps.imgName, band=inps.band,
        maskArgs=inps.maskArgs,
        verbose=inps.verbose)

    # Query profile
    prof.query_profile(inps.queryFile, inps.profWidth)

    # Plot if requested
    if inps.plot == True:
        prof.plot(cmap=inps.cmap, cbarOrient=inps.cbarOrient,
            minPct=inps.minPct, maxPct=inps.maxPct)

    # Save to file
    if inps.outName is not None:
        prof.save(inps.outName, overwrite=inps.overwrite)


    plt.show()