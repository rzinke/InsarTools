#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Apply a specified filter to a GDAL compatible image.

FUTURE IMPROVEMENTS


TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import confirm_outdir, confirm_outname_ext, load_gdal_dataset, save_gdal_dataset
from Masking import create_mask
from GeoFormatting import DS_to_extent
from Filtering import mean_filter
from Viewing import plot_raster


### PARSER ---
Description = '''Apply a specified filter to a GDAL compatible image.'''

Examples = ''''''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='imgFile', type=str,
        help='Image file name.')
    InputArgs.add_argument('-b','--band', dest='bandNb', type=int, default=1,
        help='Image band number to display.')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')
    InputArgs.add_argument('-f','--filter-type', dest='filtType', type=str, required=True,
        help='Filter type (mean).')
    InputArgs.add_argument('-w','--kernel-width', dest='kernelWidth', type=int, required=True,
        help='Kernel width (pixels).')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true',
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot results.')
    OutputArgs.add_argument('-o','--outname', dest='outName', type=str, default='Filt',
        help='GPS-corrected field.')

    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### FILTERING ---
def format_image(img, mask, verbose=False):
    '''
    Format the image before filtering to replace NaNs and masked values with
     zeros.
    '''
    if verbose == True: print('Pre-formatting image')

    # Replace NaNs with zeros
    img[np.isnan(img) == 1] = 0

    # Replace mask values
    img[mask==0] = 0

    return img


def apply_filter(img, filtType, kernelWidth=None, verbose=False):
    '''
    Apply the specified filter.
    '''
    if verbose == True: print('Applying filter')

    # Check that filter type is valid
    filterTypes = ['mean']

    if filtType not in filterTypes:
        print('Filter type {:s} not in implemented filter types:'.format(filtType))
        [print('\t{:s}'.format(filterType)) for filterType in filterTypes]
        exit()

    # Check that filter type is valid
    if filtType == 'mean':
        assert kernelWidth is not None, 'Kernel width must be specified for mean filter.'

        fimg = mean_filter(img, kernelWidth, verbose=verbose)

    return fimg



### PLOTTING ---
def plot_images(img, fimg, mask, extent, verbose=False):
    '''
    Plot the input and output images.
    '''
    if verbose == True: print('Plotting inputs and results')

    # Parameters
    cbarOrient = 'auto'

    # Plot input image
    fig0, ax0 = plot_raster(img, mask=mask, extent=extent,
        cmap='viridis', cbarOrient=cbarOrient,
        minPct=1, maxPct=99)
    ax0.set_title('Original image')

    # Plot filtered image
    fig0, ax0 = plot_raster(fimg, mask=mask, extent=extent,
        cmap='viridis', cbarOrient=cbarOrient,
        minPct=1, maxPct=99)
    ax0.set_title('Filtered image')    



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Load data
    # Load image data set
    DS = load_gdal_dataset(inps.imgFile, verbose=inps.verbose)
    extent = DS_to_extent(DS, verbose=inps.verbose)
    img = DS.GetRasterBand(inps.bandNb).ReadAsArray()

    # Create mask
    mask = create_mask(img, inps.maskArgs)

    # Pre-format image
    img = format_image(img, mask, verbose=inps.verbose)


    ## Filtering
    fimg = apply_filter(img, inps.filtType, kernelWidth=inps.kernelWidth, verbose=inps.verbose)


    ## Outputs
    # Save to file
    outName = confirm_outname_ext(inps.outName, ext=['tif'])
    confirm_outdir(inps.outName)
    save_gdal_dataset(outName, fimg, mask=mask, exDS=DS, verbose=inps.verbose)

    # Plot if requested
    if inps.plot == True: plot_images(img, fimg, mask, extent, verbose=inps.verbose)


    plt.show()