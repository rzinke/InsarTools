#!/usr/bin/env python3
'''
SHORT DESCRIPTION
View a single image.

FUTURE IMPROVEMENTS
    * Save PNG of image (without displaying)

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import load_gdal_dataset, confirm_outdir, confirm_outname_ext
from GeoFormatting import DS_to_extent
from Masking import create_mask
from Viewing import plot_raster, plot_histogram, image_clip_values, image_stats


### PARSER ---
Description = '''View a single plot of a GDAL-readable image.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='imgName',
        help='Name of file to display.')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')
    InputArgs.add_argument('-b','--band', dest='bandNb', type=int, default=1,
        help='Image band number to display.')


    DisplayArgs = parser.add_argument_group('DISPLAY PARAMS')
    DisplayArgs.add_argument('-c','--cmap', dest='cmap', type=str, default='viridis',
        help='Colormap ([viridis]).')
    DisplayArgs.add_argument('-co','--colorbar-orientation', dest='cbarOrient', type=str, default='auto',
        help='Colorbar orientation ([auto], horizontal, vertical).')
    DisplayArgs.add_argument('-minPct','--min-percent', dest='minPct', type=float, default=None,
        help='Minimum percent clip value ([None]).')
    DisplayArgs.add_argument('-maxPct','--max-percent', dest='maxPct', type=float, default=None,
        help='Maximum percent clip value ([None]).')
    DisplayArgs.add_argument('-vmin','--min-value', dest='vmin', type=float, default=None,
        help='Minimum clip value ([None]).')
    DisplayArgs.add_argument('-vmax','--max-value', dest='vmax', type=float, default=None,
        help='Maximum clip value ([None]).')
    DisplayArgs.add_argument('-eq','--equalize', dest='equalize', action='store_true',
        help='Equalize')


    HistogramArgs = parser.add_argument_group('HISTOGRAM')
    HistogramArgs.add_argument('-hist','--histogram', dest='plotHist', action='store_true')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-o','--outname', dest='outName', type=str, default=None,
        help='Output name to save image.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### LOADING AND FORMATTING ---
def load_image(imgName, bandNb=1, verbose=False):
    '''
    Load the image as a GDAL data set. Extract the specified band and geographic
     extent.
    '''
    if verbose == True: print('Loading and formatting image')

    # Load GDAL data set
    DS = load_gdal_dataset(imgName, verbose=verbose)

    # Extract image(s)
    img = DS.GetRasterBand(bandNb).ReadAsArray()

    # Replace NaN values with zeros
    img[np.isnan(img) == 1] = 0

    # Geographic extent
    extent = DS_to_extent(DS, verbose=verbose)

    return img, extent



### SAVING ---
def save_image(outName, fig, verbose=False):
    '''
    Save the image figure to an image file (PNG).
    '''
    # Confirm output directory
    confirm_outdir(outName)

    # Confirm outname extension
    outName = confirm_outname_ext(outName, ['png', 'PNG'])

    # Save figure
    fig.savefig(outName)

    # Report if requested
    if verbose == True: print('Saved to: {:s}'.format(outName))



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Load image
    # Load and extract image
    img, extent = load_image(imgName=inps.imgName,
        bandNb=inps.bandNb,
        verbose=inps.verbose)

    # Build mask
    mask = create_mask(img, inps.maskArgs, verbose=inps.verbose)


    ## Plot image
    imgFig, imgAx = plot_raster(img, mask=mask, extent=extent,
        cmap=inps.cmap, cbarOrient=inps.cbarOrient,
        vmin=inps.vmin, vmax=inps.vmax,
        minPct=inps.minPct, maxPct=inps.maxPct,
        equalize=inps.equalize)


    ## Histogram
    # Create and plot histogram if requested
    if inps.plotHist == True:
        histFig, histAx = plot_histogram(img, mask)
        histAx.set_xlim(image_clip_values(img,
            vmin=inps.vmin, vmax=inps.vmax,
            minPct=inps.minPct, maxPct=inps.maxPct))


    ## Statistics
    # Report stats if requested
    if inps.verbose == True:
        image_stats(img, mask=mask, verbose=True)


    ## Save image to file
    if inps.outName is not None:
        save_image(inps.outName, imgFig, verbose=inps.verbose)


    plt.show()