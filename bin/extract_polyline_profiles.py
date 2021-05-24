#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Extract profiles from a map that are perpendicular to a given polyline.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import os
import matplotlib.pyplot as plt
from IOsupport import load_polyline, load_gdal_dataset, confirm_outdir, confirm_outname_ext, append_fname, save_profile_data
from GeoFormatting import transform_to_extent, lola_to_xy
from Masking import create_mask
from Profiling import find_profile_geometries, extract_profile
from Viewing import plot_raster, plot_profile


### PARSER ---
Description = '''Extract profiles from a map that are perpendicular to a given polyline.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument('-l','--polyline', dest='polylineName', type=str, required=True,
        help='Polyline file with coordinates (x, y).')
    InputArgs.add_argument('-f','--image', dest='imgName', type=str, required=True,
        help='Image file.')
    InputArgs.add_argument('-b','--band', dest='band', type=int, default=1,
        help='Image band number ([1]).')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')

    ProfileArgs = parser.add_argument_group('PROFILE ARGUMENTS')
    ProfileArgs.add_argument('-w','--profile-width', dest='profWidth', type=float, default=None, required=True,
        help='Profile width in map units.')
    ProfileArgs.add_argument('-len','--profile-length', dest='profLen', type=float, default=None, required=True,
        help='Profile length in map units.')

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
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name prefix.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)


### POLYLINE PROFILE CLASS ---
class polylineProfiles:
    def __init__(self, polylineName, imgName, profWidth, profLen, band=1, maskArgs=None, verbose=False):
        '''
        Extract profiles along a polyline.
        '''
        # Parameters
        self.profWidth = profWidth
        self.profLen = profLen
        self.verbose = verbose

        # Load polyline
        self.lineX, self.lineY = load_polyline(polylineName, verbose=self.verbose)

        # Load image dataset
        self.__load_data__(imgName, band)

        # Create mask
        self.mask = create_mask(self.img, maskArgs, verbose=self.verbose)

        # Build profile geometries
        self.__build_geometries__()

        # Extract profile
        self.__extract_profiles__()


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
    def __build_geometries__(self):
        '''
        Find profile geometries.
        Return a list of profile_geometry objects, each containing the geometry
         information of a single profile.
        '''
        # Build geometries
        self.profGeoms = find_profile_geometries(lineX=self.lineX, lineY=self.lineY,
            profWidth=self.profWidth, profLen=self.profLen,
            verbose=self.verbose)

        # Number of profiles
        self.nProfiles = len(self.profGeoms)


    def __extract_profiles__(self):
        '''
        Convert units to pixel coordinates and extract profiles.
        '''
        # Setup
        self.profDists = []
        self.profPts = []

        for profGeom in self.profGeoms:
            # Profile start and end points
            xStart = np.array(profGeom.profStart[0])
            yStart = np.array(profGeom.profStart[1])

            xEnd = np.array(profGeom.profEnd[0])
            yEnd = np.array(profGeom.profEnd[1])

            # Convert to pixel coordinates
            pxStart, pyStart = lola_to_xy(self.tnsf, xStart, yStart)
            pxEnd, pyEnd = lola_to_xy(self.tnsf, xEnd, yEnd)

            # Extract profiles
            pxDist, profPts = extract_profile(self.img,
                pxStart=pxStart[0], pyStart=pyStart[0],
                pxEnd=pxEnd[0], pyEnd=pyEnd[0],
                pxWidth=self.profWidth/self.tnsf[1],
                mask=self.mask,
                verbose=self.verbose)

            # Scale profile distance from pixels to map units
            pxSize = self.tnsf[1]
            profDist = pxDist*pxSize

            # Append to lists
            self.profDists.append(profDist)
            self.profPts.append(profPts)


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

        # Plot polyline
        axMap.plot(self.lineX, self.lineY, 'k')

        # Plot profile anchor points
        for profGeom in self.profGeoms:
            # Plot mid points
            axMap.plot(profGeom.midAnchor[0], profGeom.midAnchor[1], 'k+')

            # Plot profile corners
            plot_profile(profGeom, fig=mapFig, ax=axMap)

        # Plot profiles
        profFig, profAx = plt.subplots(figsize=(9,5))

        for i in range(self.nProfiles):
            profAx.plot(self.profDists[i], self.profPts[i],
                marker='.', linewidth=0, label='Profile {:d}'.format(i+1))

        profAx.legend()


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
        for i in range(self.nProfiles):
            # Append profile name
            fname = append_fname(outName, '_{:d}'.format(i+1))

            # Report if requested
            if self.verbose == True: print('... {:s}'.format(os.path.basename(fname)))

            # Save to file
            save_profile_data(fname, self.profGeoms[i].profStart, self.profGeoms[i].profEnd,
                self.profDists[i], self.profPts[i])



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Compute profiles
    profs = polylineProfiles(polylineName=inps.polylineName, imgName=inps.imgName,
        maskArgs=inps.maskArgs,
        profWidth=inps.profWidth, profLen=inps.profLen,
        verbose=inps.verbose)

    # Plot if requested
    if inps.plot == True:
        profs.plot(cmap=inps.cmap, cbarOrient=inps.cbarOrient,
            minPct=inps.minPct, maxPct=inps.maxPct)

    # Save to file
    if inps.outName is not None:
        profs.save(inps.outName)


    plt.show()