#!/usr/bin/env python3
'''
SHORT DESCRIPTION
View Relax displacement outputs.

FUTURE IMPROVEMENTS
    * Zoom in on fault

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import detect_relax_files, load_gdal_datasets
from GeoFormatting import DS_to_extent
from Masking import create_mask
from Viewing import dataset_clip_values, raster_multiplot
from project_to_los import LOSproject


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
    DisplayArgs.add_argument('-vmin','--min-value', dest='vmin', type=float, default=None,
        help='Minimum clip value ([None]).')
    DisplayArgs.add_argument('-vmax','--max-value', dest='vmax', type=float, default=None,
        help='Maximum clip value ([None]).')


    ProjectionArgs = parser.add_argument_group('PROJECTION')
    ProjectionArgs.add_argument('--projection-convention', dest='projConvention', type=str, default=None,
        help='Projection convention ([None], ARIA, ISCE].')
    ProjectionArgs.add_argument('--incidence', dest='incInpt', type=str, default=None,
        help='Incidence (float or filename).')
    ProjectionArgs.add_argument('--azimuth', dest='azInpt', type=str, default=None,
        help='Azimuth (float or filename).')
    ProjectionArgs.add_argument('--geometry', dest='geomFile', type=str, default=None,
        help='ISCE geometry file.')
    ProjectionArgs.add_argument('--wrap', dest='wrap', type=float, default=None,
        help='Wrap signal to value ([None]).')


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

    # Number of epochs
    nEdata = len(Edata)
    nNdata = len(Ndata)
    nUdata = len(Udata)
    assert nEdata == nNdata == nUdata, \
        'Number of E ({:d}), N ({:d}), U ({:d}) data are not the same'.format(nEdata, nNdata, nUdata)
    nEpochs = nEdata

    # Extract images from data sets
    Eimgs = [DS.GetRasterBand(1).ReadAsArray() for DS in Edata.values()]
    Nimgs = [DS.GetRasterBand(1).ReadAsArray() for DS in Ndata.values()]
    Uimgs = [DS.GetRasterBand(1).ReadAsArray() for DS in Udata.values()]

    # Determine clip values
    Emin, Emax = dataset_clip_values(Eimgs, verbose=inps.verbose)
    Nmin, Nmax = dataset_clip_values(Nimgs, verbose=inps.verbose)
    Umin, Umax = dataset_clip_values(Uimgs, verbose=inps.verbose)

    if inps.vmin is not None:
        vmin = inps.vmin
    else:
        vmin = np.min([Emin, Nmin, Umin])

    if inps.vmax is not None:
        vmax = inps.vmax
    else:
        vmax = np.max([Emax, Nmax, Umax])

    # Geographic extent
    extent = DS_to_extent(list(Edata.values())[0])  # same geographic extent


    ## Plotting
    # Loop through epochs
    for i in range(nEpochs):
        # Retrieve images for the 
        Eimg = Eimgs[i]
        Nimg = Nimgs[i]
        Uimg = Uimgs[i]

        # Format inputs based on projection or not
        if inps.projConvention is not None:
            # Project into LOS
            Pimg = LOSproject(convention=inps.projConvention,
                incInpt=inps.incInpt, azInpt=inps.azInpt, geomFile=inps.geomFile,
                Einpt=Eimg, Ninpt=Nimg, Vinpt=Uimg,
                verbose=inps.verbose)

            LOS = Pimg.LOS

            # Wrap if requested
            if inps.wrap is not None:
                LOS %= inps.wrap
                LOS *= vmax/inps.wrap

            # Format images as a stack
            imgs = [Eimg, Nimg, Uimg, LOS]

            # Format titles
            titles = ['E', 'N', 'U', 'Proj']

        else:
            # Format images as a stack
            imgs = [Eimg, Nimg, Uimg]

            # Format titles
            titles = ['E', 'N', 'U']

        # Format super title
        suptitle = Enames[i][:3]

        # Spawn figure
        fig, axes = plt.subplots(figsize=(8,6), nrows=2, ncols=2)

        # Plot model results
        fig, axes = raster_multiplot(imgs, mrows=2, ncols=2, extent=extent,
            cmap=inps.cmap, cbarOrient=inps.cbarOrient,
            vmin=vmin, vmax=vmax,
            titles=titles, suptitle=suptitle,
            fig=fig, axes=axes)


    plt.show()