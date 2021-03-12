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
import matplotlib.pyplot as plt
from IOsupport import load_gdal_dataset, confirm_outdir, confirm_overwrite
from GeoFormatting import transform_to_extent, lola_to_xy
from Masking import create_mask
from Profiling import extract_profile
from Viewing import plot_raster


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
    ProfileArgs.add_argument('-q','--queryLoLa', dest='qLoLa', nargs=2, required=True,
        help='Query point in geographic coordinates (e.g., 89.25,34.07 88.72,37.02)')
    ProfileArgs.add_argument('-w','--profile-width', dest='profWidth', type=float, default=None,
        help='Profile width in pixels.')

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

        # Extract geographic information
        self.M, self.N = DS.RasterYSize, DS.RasterXSize
        self.tnsf = DS.GetGeoTransform()
        self.extent = transform_to_extent(self.tnsf, self.M, self.N, verbose=self.verbose)        


    ## Profiling
    def query_profile(self, qLoLa, width):
        '''
        Query the profile between the points given in qLoLa and with the
         specified width.
        '''
        if self.verbose == True: print('Querying location:')

        # Parse profile points
        self.__parse_prof_points__(qLoLa)

        # Format profile width
        self.__format_width__(width)

        # Extract profile
        self.profDist, self.profPts = extract_profile(img=self.img,
            pxStart=self.pxStart, pyStart=self.pyStart,
            pxEnd=self.pxEnd, pyEnd=self.pyEnd,
            width=self.width,
            mask=self.mask,
            verbose=self.verbose)

        # Number of profile points
        self.nPts = len(self.profPts)


    def __parse_prof_points__(self, qLoLa):
        '''
        Confirm that the query points are in the correct format, i.e.,
            qLoLa = ((startLon, startLat) ())
        and are within the extent of the map.
        '''
        # Format query points
        assert len(qLoLa) == 2, 'Not enough query points provided'
        for i in range(2):
            assert len(qLoLa[i].split(',')) == 2, 'Query point must have format (x, y)'

        # Parse lon/lat points
        startLoLa, endLoLa = qLoLa

        startLon, startLat = startLoLa.split(',')
        endLon, endLat = endLoLa.split(',')

        self.startLon = float(startLon)
        self.startLat = float(startLat)

        self.endLon = float(endLon)
        self.endLat = float(endLat)


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


    def __format_width__(self, width):
        '''
        Format the profile with based on the pixel size if no set width is
         specified.
        '''
        # If width is not specified, use pixel size
        if width is None:
            self.width = 2  # use x-step
        else:
            self.width = width/self.tnsf[1]  # use specfied width in map units


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

        # Plot profile start
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
        if outName[-4:] != '.txt': outName += '.txt'

        # File formatting
        metadata = 'start: {:f} {:f}\nend: {:f} {:f}\n'
        header = '# distance amplitude\n'
        dataStr = '{:f} {:f}\n'

        # Check overwrite status
        if overwrite == False:
            overwrite = confirm_overwrite(outName)

        if overwrite == True:
            with open(outName, 'w') as profFile:
                profFile.write(metadata)
                profFile.write(header)
                for i in range(self.nPts):
                    profFile.write(dataStr.format(self.profDist[i], self.profPts[i]))
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
    prof.query_profile(inps.qLoLa, inps.profWidth)

    # Plot if requested
    if inps.plot == True:
        prof.plot(cmap=inps.cmap, cbarOrient=inps.cbarOrient,
            minPct=inps.minPct, maxPct=inps.maxPct)

    # Save to file
    if inps.outName is not None:
        prof.save(inps.outName, overwrite=inps.overwrite)


    plt.show()