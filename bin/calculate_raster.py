#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Create a new raster by doing math on existing rasters.

INHERITANCES
IOsupport: load_gdal_dataset, save_gdal_dataset, pick_dataset
Masking: create_mask
RasterResampling: match_rasters
GeoResampling: DS_to_extent
Viewing: image_percentiles

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import load_gdal_datasets, pick_dataset, check_outname_ext, save_gdal_dataset
from Masking import create_mask
from RasterResampling import match_rasters
from GeoFormatting import DS_to_extent
from Viewing import image_percentiles


### PARSER ---
Description = '''Enter a formula to modify and save a raster map. This function works on multiple data sets.'''

Examples = ''''''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    # Necessary 
    parser.add_argument(dest='varStr', nargs='+', type=str,
        help='Image designation/filename pairs, e.g., A img1.tif B img2.tif C img3.tif ...')
    parser.add_argument('-f','--formula', dest='formula', type=str, required=True,
        help='Math operation(s) to be applied. Refer to the images as \"A\"-\"Z\".')
    parser.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values.')

    # Outputs
    parser.add_argument('-v','--verbose', dest='verbose', action='store_true',
        help='Verbose mode.')
    parser.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot outputs.')
    parser.add_argument('--plot-inputs', dest='plotInputs', action='store_true',
        help='Plot inputs.')
    parser.add_argument('-o','--outName', dest='outName', default='Out',
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



### MASKING ---
def build_mask(datasets, maskArgs, verbose=False):
    '''
    Build a master mask from the given datasets and mask arguments.
    Data sets must be resampled to the same dimensions.
    '''
    if verbose == True: print('Creating common mask')

    # Setup
    varNames = list(datasets.keys())  # variable names
    M = datasets[varNames[0]].RasterYSize  # map y dim
    N = datasets[varNames[0]].RasterXSize  # map x dim
    commonMask = np.ones((M, N))  # empty common mask

    # Confirm that map sizes are the same
    for varName in varNames:
        # Retrieve image values
        img = datasets[varName].GetRasterBand(1).ReadAsArray()

        # Check that img size is the same for each data set
        assert img.shape == (M, N), 'Dataset {:s} is not the same size as those before it'

        # Compute individual mask
        mask = create_mask(img, maskArgs, verbose=verbose)

        # Modify common mask
        commonMask[mask==0] = 0

    return commonMask



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

    # Loop through inputs to plot
    for i, varName in enumerate(datasets.keys()):
        # Grab extent from first data set
        if i == 0: extent = DS_to_extent(datasets[varName])

        # Retrieve image
        img = datasets[varName].GetRasterBand(1).ReadAsArray()

        # Mask image
        img = np.ma.array(img, mask=(mask==0))

        # Compute clip values
        vmin, vmax = image_percentiles(img)

        # Plot image
        cax = axes[i].imshow(img, extent=extent,
            cmap='jet', vmin=vmin, vmax=vmax)

        # Format plot
        axes[i].set_title(varName)
        fig.colorbar(cax, ax=axes[i], orientation='horizontal')

    # Plot common mask
    cax = axes[-1].imshow(mask, extent=extent)

    # Format mask plot
    axes[-1].set_title('Mask')
    fig.colorbar(cax, ax=axes[-1], orientation='horizontal')
    fig.suptitle('Inputs')


def plot_output(img, mask, extent):
    '''
    Plot the output dataset.
    '''
    # Spawn figure
    fig, ax = plt.subplots()

    # Mask image
    img = np.ma.array(img, mask=(mask==0))

    # Compute clip values
    vmin, vmax = image_percentiles(img)

    # Plot image
    cax = ax.imshow(img, extent=extent,
        cmap='jet', vmin=vmin, vmax=vmax)

    # Format plot
    ax.set_title('Result')
    fig.colorbar(cax, ax=ax, orientation='horizontal')



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
    mask = build_mask(datasets, inps.maskArgs, verbose=inps.verbose)

    # Plot inputs if requested
    if inps.plotInputs == True:
        plot_inputs(datasets, mask)


    ## Evaluate formula
    if inps.verbose == True: print('*'*32)  # stdout break

    # Compute image
    img = apply_formula(datasets, inps.formula, mask=None, verbose=inps.verbose)


    ## Outputs
    # Save as GDAL dataset
    outName = check_outname_ext(inps.outName)
    save_gdal_dataset(outName, img, pick_dataset(datasets), verbose=inps.verbose)

    # Plot output
    if inps.plot == True:
        plot_output(img, mask, DS_to_extent(pick_dataset(datasets)))


    plt.show()