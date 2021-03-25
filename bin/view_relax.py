#!/usr/bin/env python3
'''
SHORT DESCRIPTION
View Relax displacement outputs.

FUTURE IMPROVEMENTS
    * Save PNGs of images (without displaying)

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import matplotlib.pyplot as plt
from IOsupport import detect_relax_files, load_gdal_datasets
from GeoFormatting import DS_to_extent
from Masking import create_mask
from Viewing import raster_multiplot


### PARSER ---
Description = '''View Relax displacement outputs.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='relaxDir', type=str,
        help='Relax file directory.')
    InputArgs.add_argument('-b','--band', dest='bandNb', type=int, default=1,
        help='Image band number to display.')


    DisplayArgs = parser.add_argument_group('DISPLAY PARAMS')
    DisplayArgs.add_argument('-c','--cmap', dest='cmap', type=str, default='viridis',
        help='Colormap ([viridis]).')
    DisplayArgs.add_argument('-co','--colorbar-orientation', dest='cbarOrient', type=str, default='vertical',
        help='Colorbar orientation (horizontal, [vertical]).')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-o','--outname', dest='outName', type=str, default=None,
        help='Output name to save image.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Load data
    # Detect displacement filenames
    Enames, Nnames, Unames = detect_relax_files(inps.relaxDir, verbose=inps.verbose)

    # Load displacement data sets
    Edata = load_gdal_datasets(Enames, verbose=inps.verbose)
    Ndata = load_gdal_datasets(Nnames, verbose=inps.verbose)
    Udata = load_gdal_datasets(Unames, verbose=inps.verbose)

    # List keys
    Ekeys = list(Edata.keys())
    Nkeys = list(Ndata.keys())
    Ukeys = list(Udata.keys())

    # Check that data sets are equal length
    assert len(Ekeys) == len(Nkeys) == len(Ukeys), \
        'E, N, and up data sets must be equal length'

    nDatasets = len(Ekeys)


    ## Plot data
    for i in range(nDatasets):
        # Format images into list
        imgs = [Edata[Ekeys[i]], Ndata[Nkeys[i]], Udata[Ukeys[i]]]

        # Plot images
        fig, axes = raster_multiplot(imgs, mrows=2, ncols=2,
            cmap=inps.cmap, cbarOrient=inps.cbarOrient)

        axes[0].set_title('East')
        axes[1].set_title('North')
        axes[2].set_title('Up')


    plt.show()