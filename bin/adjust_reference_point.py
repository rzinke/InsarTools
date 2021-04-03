#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Change the reference point of a raster data set.

FUTUTRE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import load_gdal_dataset, confirm_outdir, confirm_outname_ext, save_gdal_dataset
from Masking import mask_dataset
from GeoFormatting import DS_to_extent, lola_to_xy
from Checks import check_loc_in_dataset
from Viewing import image_percentiles


### PARSER ---
Description = '''Change the reference point of a raster data set.'''

Examples = ''''''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    # Inputs
    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='imgName', type=str,
        help='Image name.')
    InputArgs.add_argument('-lon','--longitude', dest='lon', type=float, default=None,
        help='Reference longitude.')
    InputArgs.add_argument('-lat','--latitude', dest='lat', type=float, default=None,
        help='Reference latitude.')
    InputArgs.add_argument('-r','--radius', dest='radius', type=int, default=0,
        help='\"Search radius\" of pixels around the given reference pixel. ([0], 1, 2, ...)')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')
    InputArgs.add_argument('-i','--interactive', dest='interactive', action='store_true',
        help='Interactive mode.')

    # Outputs
    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true',
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot outputs.')
    OutputArgs.add_argument('-o','--outname', dest='outName', type=str, default='Out',
        help='Output head name.')

    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### REFERENCE POINT ---
class adjust_ref_pt:
    ## Defaults
    interactive = False  # declare static mode as default


    ## Setup
    def __init__(self, imgName, maskArgs, verbose=False):
        '''
        Adjust the reference point of a given data set to the new specified lon and lat.
        First, convert the lon/lat coords to pixel values and check that they are not masked.
        Then, adjust the velocity field to the new reference point by subtracting the velocity
         value at that location.

        Initialize by loading the original data set and creating a mask.
        '''
        # Parameters
        self.verbose = verbose

        # Load data set
        if self.verbose == True: print('*'*32)
        self.DS = load_gdal_dataset(imgName, verbose=self.verbose)

        # Determine geographic information
        self.tnsf = self.DS.GetGeoTransform()
        self.extent = DS_to_extent(self.DS)

        # Create mask
        self.mask = mask_dataset(self.DS, maskArgs, verbose=self.verbose)

        # Retrieve image
        self.__retrieve_image__()


    def __retrieve_image__(self):
        '''
        Retrieve original image from GDAL data set, mask, and copy to referenced image.
        '''
        # Retrieve original image
        self.img = self.DS.GetRasterBand(1).ReadAsArray()

        # Mask image
        self.img = np.ma.array(self.img, mask=(self.mask==0))

        # Copy to re-referencing image
        self.refImg = self.img.copy()


    ## Adjustment
    def __adjust_to_pixel__(self, px, py, radius):
        '''
        Adjust the image to the given pixel location.
        '''
        # Determine value at reference point
        self.__get_ref_value__(px, py, radius)

        # Remove reference point from data set
        self.refImg = self.img - self.refValue

        # Report if specified
        if self.verbose == True:
            print('Removed reference value: {:f}'.format(self.refValue))


    def __get_ref_value__(self, px, py, radius):
        '''
        Get the reference value from the original image.
        '''
        # Range of pixels based on location and radius
        yRange = np.arange(py-radius, py+radius+1)
        xRange = np.arange(px-radius, px+radius+1)
        
        # Adjustment value
        refValues = self.img[yRange,xRange]
        self.refValue = np.median(refValues.compressed().flatten())


    ## Plotting
    def __plot_adjusted_image__(self, fig, ax, vmin, vmax):
        '''
        Plot the re-referenced image.
        '''
        # Plot image
        cax = ax.imshow(self.refImg, extent=self.extent,
            cmap='jet', vmin=vmin, vmax=vmax)

        return cax


    def __plot_ref_point__(self, fig, ax, lon, lat):
        '''
        Plot the current reference point.
        '''
        self.ax.scatter(lon, lat, facecolor='w', edgecolor='k')


    ## Standard mode
    def adjust_reference(self, lon, lat, radius=0):
        '''
        Adjust the image to the new reference point using non-GUI adjustment.
        '''
        # Coordinate values
        self.lon = lon
        self.lat = lat

        # Check that reference point is within valid pixels; convert lon/lat to pixel coordinates
        px, py = check_loc_in_dataset(self.DS, self.mask, self.lon, self.lat, verbose=self.verbose)

        # Adjust image
        self.__adjust_to_pixel__(px, py, radius)


    def plot(self):
        '''
        Plot results of reference point adjustment from the standard workflow.
        '''
        # Spawn figure
        fig, ax = plt.subplots()
        cbarOrient = 'vertical'

        # Determine clip values
        vmin, vmax = image_percentiles(self.refImg)

        # Plot re-referenced image
        cax = self.__plot_adjusted_image__(fig, ax, vmin, vmax)

        # Plot reference point
        ax.scatter(self.lon, self.lat, facecolor='w', edgecolor='k')

        # Format plot
        fig.colorbar(cax, ax=ax, orientation=cbarOrient)


    ## Interactive mode
    def plot_interactive(self, lon=None, lat=None, radius=0):
        '''
        Interactive plot for reference point adjustment.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Interactive reference point adjustment')

        # Record search radius
        self.radius = radius

        # Find clip values
        self.vmin, self.vmax = image_percentiles(self.img)

        # Create interactive plot
        self.fig, self.ax = plt.subplots()
        cbarOrient = 'vertical'

        # Adjust to pre-specified reference point
        if lon is not None and lat is not None:
            self.adjust_reference(lon, lat, radius)
            self.__plot_ref_point__(self.fig, self.ax, lon, lat)

        # Plot original image
        cax = self.__plot_adjusted_image__(self.fig, self.ax, self.vmin, self.vmax)

        # Colorbar
        self.fig.colorbar(cax, ax=self.ax, orientation=cbarOrient)

        # Interact with image
        self.fig.canvas.mpl_connect('button_press_event', self.__live_adjust__)


    def __live_adjust__(self, event):
        '''
        Perform a "live" reference point adjustment.
        '''
        # Record click values
        lon = event.xdata
        lat = event.ydata

        # Convert to pixel values
        px, py = lola_to_xy(self.tnsf, lon, lat, verbose=self.verbose)

        # Check that mask is not selected
        if self.mask[py,px] == 1:
            # ... Continue with valid value

            # Adjust to selected pixel
            self.__adjust_to_pixel__(px, py, self.radius)

            # Clear axis
            self.ax.cla()

            # Replot image
            self.__plot_adjusted_image__(self.fig, self.ax, self.vmin, self.vmax)

            # Plot reference point
            self.__plot_ref_point__(self.fig, self.ax, lon, lat)

            # Update image
            self.fig.canvas.draw()

        else:
            # ... Report invalid value
            print('Location ({:f} {:f}) is masked, try again'.format(lon, lat))



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Adjust reference point
    # Instantiate reference adjustment class
    refAdjust = adjust_ref_pt(inps.imgName, inps.maskArgs, verbose=inps.verbose)

    if inps.interactive == True:
        # Interactive mode
        refAdjust.plot_interactive(lon=inps.lon, lat=inps.lat)
        plt.show()

    else:
        # Standard mode
        refAdjust.adjust_reference(lon=inps.lon, lat=inps.lat, radius=inps.radius)

        if inps.plot == True:
            refAdjust.plot()
            plt.show()


    ## Save data set
    # Checks
    outName = confirm_outname_ext(inps.outName, ['tif'])  # confirm output extension
    confirm_outdir(outName)  # confirm output directory exists

    # Save data set
    save_gdal_dataset(outName, refAdjust.refImg,
        mask=refAdjust.mask, exDS=refAdjust.DS,
        verbose=inps.verbose)