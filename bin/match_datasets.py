#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Resample several rasters to a common extent. Essentially a wrapper for match_rasters.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import os
import argparse
import matplotlib.pyplot as plt
from IOsupport import confirm_outdir, rel_to_abs_paths, load_gdal_datasets, append_fname, save_gdal_dataset
from RasterResampling import match_rasters
from Checks import check_dataset_sizes
from Masking import mask_datasets
from Viewing import plot_raster


### PARSER ---
Description = '''Resample several rasters to a common geographic extent.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='imgFiles', type=str, nargs='+',
        help='Files to be resampled.')
    InputArgs.add_argument('-c','--cropping', dest='cropping', type=str, default='union',
        help='Cropping scheme ([union], intersection).')
    InputArgs.add_argument('-r','--resolution', dest='resolution', type=str, default='fine',
        help='Resolution ([fine], coarse).')

    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name prefix. ([Out]).')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true', 
        help='Plot inputs and outputs.')
    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### NAME FORMATTING ---
def format_outnames(inNames, prefix, verbose=False):
    '''
    Format the output names by adding the output name to the front, and "resampled" to the back.
    '''
    # Setup
    outNames = []

    # Loop through input names
    for inName in inNames:
        # Copy input name
        outName = inName[:]

        # Add output to front
        basename = os.path.basename(outName)
        basename = '{:s}_{:s}'.format(prefix, basename)
        dirname = os.getcwd()
        outName = os.path.join(dirname, basename)

        # Add "resampled" appendix to the back
        outName = append_fname(outName, '_res', verbose=True)

        # Append to list
        outNames.append(outName)

    # Confirm output directory
    confirm_outdir(outNames[0])

    return outNames



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()

    # Convert relative paths to absolute
    imgFiles = rel_to_abs_paths(inps.imgFiles)

    # Load data sets
    datasets = load_gdal_datasets(imgFiles, verbose=inps.verbose)


    ## Resample
    datasets = match_rasters(datasets, cropping=inps.cropping, resolution=inps.resolution,
        verbose=inps.verbose)


    ## Save data
    # Format output names
    outNames = format_outnames(imgFiles, inps.outName, verbose=inps.verbose)

    # Save each resampled data set
    for i, dsName in enumerate(datasets.keys()):
        # Isolate data set
        DS = datasets[dsName]

        # Save data set
        save_gdal_dataset(outNames[i], DS.GetRasterBand(1).ReadAsArray(), exDS=DS, verbose=inps.verbose)


    ## Plotting
    if inps.plot == True:
        # Mask
        masks = mask_datasets(datasets, ['bg'])

        # Plot
        for i, dsName in enumerate(datasets.keys()):
            fig, ax = plot_raster(datasets[dsName], mask=masks)
            ax.set_title(os.path.basename(dsName))


    plt.show()