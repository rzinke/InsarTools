#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Create a new raster by doing math on existing rasters.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import load_gdal_datasets, pick_dataset, confirm_outdir, confirm_outname_ext, save_gdal_dataset
from Masking import create_common_mask
from RasterResampling import match_rasters
from GeoFormatting import DS_to_extent
from Viewing import plot_raster


### PARSER ---
Description = '''Enter a formula to modify and save a raster map. This function works on multiple data sets.'''

Examples = ''''''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    # Inputs
    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='varStr', nargs='+', type=str,
        help='Image designation/filename pairs, e.g., A img1.tif B img2.tif C img3.tif ...')
    InputArgs.add_argument('-f','--formula', dest='formula', type=str, required=True,
        help='Math operation(s) to be applied. Refer to the images as \"A\"-\"Z\".')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')

    # Outputs
    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true',
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot outputs.')
    OutputArgs.add_argument('--plot-inputs', dest='plotInputs', action='store_true',
        help='Plot inputs.')
    OutputArgs.add_argument('-o','--outName', dest='outName', default='Out',
        help='Name of output file.')

    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### LOADING ---
def load_images(varStr, verbose=False):
    '''
    Load images into a dictionary, such that the dict key is the variable name, and the value is the map.
    '''
    if verbose == True: print('Loading data sets')

    # Parse variable names from file names
    varNames = varStr[0::2]  # variable names
    fileNames = varStr[1::2]  # file names

    # Load files
    datasets = load_gdal_datasets(fileNames, varNames, verbose=verbose)

    return datasets



### OPERATIONS ---
def apply_formula(datasets, formula, mask, verbose=False):
    '''
    Format and mask images, then apply the specified formula.
    '''
    if verbose == True: print('Applying formula ...')

    # Create dictionary of images
    imgs = {}  # empty image dictionary
    for varName in datasets.keys():
        imgs[varName] = datasets[varName].GetRasterBand(1).ReadAsArray()

    # Mask images
    for varName in imgs.keys():
        imgs[varName] = np.ma.array(imgs[varName], mask=(mask==0))

    # Apply formula
    img = eval(formula, imgs)

    # Unmask and fill nodata values
    img.mask = np.ma.nomask
    img[mask==0] = 0

    return img



### PLOTTING ---
def plot_inputs(datasets, mask):
    '''
    Plot input data sets and mask.
    '''
    # Setup
    nDatasets = len(list(datasets.keys()))
    fig, axes = plt.subplots(ncols=nDatasets+1)
    cbarOrient = 'horizontal'

    # Loop through inputs to plot
    for i, varName in enumerate(datasets.keys()):
        # Plot data set
        fig, axes[i] = plot_raster(datasets[varName], mask=mask,
            cmap='jet', cbarOrient=cbarOrient, minPct=1, maxPct=99,
            fig=fig, ax=axes[i])
        axes[i].set_title(varName)

    # Plot common mask
    fig, axes[-1] = plot_raster(mask,
        cmap='magma', cbarOrient=cbarOrient,
        vmin=0, vmax=1,
        fig=fig, ax=axes[-1])
    axes[-1].set_title('Mask')
    axes[-1].set_xticks([]); axes[-1].set_yticks([])
    fig.suptitle('Inputs')


def plot_output(img, mask, extent):
    '''
    Plot the output dataset.
    '''
    # Plot raster
    fig, ax = plot_raster(img, mask=mask, extent=extent,
        cmap='jet', cbarOrient='vertical', minPct=1, maxPct=99)

    # Format plot
    ax.set_title('Result')



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Load and format data
    if inps.verbose == True: print('*'*32)  # stdout break

    # Load data sets
    datasets = load_images(inps.varStr, verbose=inps.verbose)

    # Resample to same size
    datasets = match_rasters(datasets, cropping='intersection', verbose=inps.verbose)


    ## Determine mask
    # Create master mask for data set
    mask = create_common_mask(datasets, inps.maskArgs, verbose=inps.verbose)

    # Plot inputs if requested
    if inps.plotInputs == True:
        plot_inputs(datasets, mask)


    ## Evaluate formula
    if inps.verbose == True: print('*'*32)  # stdout break

    # Compute image
    img = apply_formula(datasets, inps.formula, mask=None, verbose=inps.verbose)


    ## Save to file
    # Checks
    confirm_outdir(inps.outName)  # confirm output directory exists
    outName = confirm_outname_ext(inps.outName, ['tif', 'tiff'])  # confirm file extension

    # Save data set
    save_gdal_dataset(outName, img, mask=mask, exDS=pick_dataset(datasets), verbose=inps.verbose)


    ## Plot
    if inps.plot == True: plot_output(img, mask, DS_to_extent(pick_dataset(datasets)))


    plt.show()