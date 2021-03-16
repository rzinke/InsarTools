#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Fit and remove an nth-degree ramp from the given data set.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import load_gdal_dataset, confirm_outname_ext, save_gdal_dataset
from Masking import create_mask
from Fitting import fit_surface
from GeoFormatting import DS_to_extent
from Viewing import image_percentiles


### PARSER ---
Description = '''Fit and remove an nth-degree ramp from the given data set.'''

Examples = ''''''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='imgFile', type=str,
        help='Image file name.')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')
    InputArgs.add_argument('-d','--degree', dest='degree', type=int, default=1,
        help='Fit an nth-degree surface. (0, [1], 2, 3, ...).')
    InputArgs.add_argument('-df','--decimation-factor', dest='decimation', type=int, default=0,
        help='Data set decimation factor, 10^decimation. ([0], 1, 2, ...).')
    InputArgs.add_argument('--keep-offset', dest='keepOffset', action='store_true',
        help='Retain offset value (z-position).')

    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true',
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot results.')
    OutputArgs.add_argument('-o','--outname', dest='outName', type=str, default='Out',
        help='Output head.')

    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### PLOTTING ---
def plot_results(OGimg, ramp, derampImg, mask, extent):
    '''
    Plot the original image, ramp, and deramped image.
    '''
    # Setup
    fig, [axOG, axRamp, axDeramp] = plt.subplots(ncols=3)
    cbarOrient = 'horizontal'

    # Plot original image
    OGimg = np.ma.array(OGimg, mask=(mask==0))  # mask image
    vmin, vmax = image_percentiles(OGimg)  # value extremes
    caxOG = axOG.imshow(OGimg, extent=extent,
        cmap='jet', vmin=vmin, vmax=vmax)

    # Format original image
    axOG.set_title('Orig. image')
    fig.colorbar(caxOG, ax=axOG, orientation=cbarOrient)

    # Plot ramp
    caxRamp = axRamp.imshow(ramp, extent=extent,
        cmap='jet', vmin=vmin, vmax=vmax)

    # Format ramp
    axRamp.set_title('Ramp')
    fig.colorbar(caxRamp, ax=axRamp, orientation=cbarOrient)

    # Plot deramped image
    derampImg = np.ma.array(derampImg, mask=(mask==0))
    vmean = np.mean([vmin, vmax])
    vmin += vmean; vmax += vmax
    caxDeramp = axDeramp.imshow(derampImg, extent=extent,
        cmap='jet', vmin=vmin, vmax=vmax)

    # Format deramped image
    axDeramp.set_title('Deramped')
    fig.colorbar(caxDeramp, ax=axDeramp, orientation=cbarOrient)



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Load data
    if inps.verbose == True:
        print('*'*32)
        print('Loading data')

    # Load data set
    DS = load_gdal_dataset(inps.imgFile, verbose=inps.verbose)

    # Image and metadata
    img = DS.GetRasterBand(1).ReadAsArray()
    M, N = img.shape

    _, dx, _, _, _, dy = DS.GetGeoTransform()
    dy = abs(dy)

    # Create mask
    mask = create_mask(img, inps.maskArgs, verbose=inps.verbose)


    ## Fit ramp
    if inps.verbose == True: print('*'*32)

    # Fit surface
    ramp, B = fit_surface(img, mask, inps.degree, dx=dx, dy=dy, decimation=inps.decimation,
        verbose=inps.verbose)


    ## Remove ramp
    # Deramp image
    derampImg = img - ramp

    # Retain offset if specified
    if inps.keepOffset == True: derampImg += B[0]

    # Reapply mask
    derampImg[mask==0] = 0


    ## Save to file
    outName = confirm_outname_ext(inps.outName, ['tif', 'tiff'], verbose=inps.verbose)
    save_gdal_dataset(outName, derampImg, mask=mask, exDS=DS, verbose=inps.verbose)


    ## Plot results
    # Plot if requested
    if inps.plot == True:
        plot_results(img, ramp, derampImg, mask, DS_to_extent(DS))


    plt.show()