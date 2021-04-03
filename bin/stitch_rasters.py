#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Stitch maps belonging to two or more spatial data sets.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import confirm_outdir, confirm_outname_ext, load_gdal_datasets, pick_dataset, save_gdal_dataset
from RasterResampling import match_rasters
from Masking import mask_datasets, find_mask_overlap
from GeoFormatting import DS_to_extent
from Checks import check_dataset_sizes
from Viewing import image_percentiles, plot_raster
from Fitting import fit_surface


### PARSER ---
Description = '''Stitch multiple georeferenced data sets.'''

Examples = ''''''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='imgFiles', nargs='+', type=str,
        help='Image file names.')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')
    InputArgs.add_argument('-a','--adjust', dest='adjust', action='store_true',
        help='Adjust subsequent maps to the first map.')
    InputArgs.add_argument('-e','--expectation', dest='expectation', type=str, default='median',
        help='Expectation operator ([\'mean\'], \'median\', \'bridge\')')

    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true',
        help='Verbose mode.')
    OutputArgs.add_argument('--save-resampled', dest='saveResampled', action='store_true',
        help='Save resampled maps.')
    OutputArgs.add_argument('--plot-inputs', dest='plotInputs', action='store_true',
        help='Plot inputs.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot results.')
    OutputArgs.add_argument('-o','--outname', dest='outName', type=str, default='Out',
        help='Output head.')

    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### STITCHING ---
def stitch_datasets(datasets, masks, adjust=False, expectation='median', verbose=False):
    '''
    Stitch two or more data sets.
    '''
    if verbose == True: print('Stitching maps')

    # Setup
    M, N = check_dataset_sizes(datasets)
    stitchedImg = np.zeros((M, N))
    stitchedMask = np.zeros((M, N))

    # Loop through data sets
    for i, dsName in enumerate(datasets.keys()):
        # Retrieve image
        img = datasets[dsName].GetRasterBand(1).ReadAsArray()

        # Retrieve mask
        mask = masks[dsName]

        # Apply mask
        img[mask==0] = 0

        # Find area of overlap between masks
        overlap = find_mask_overlap(stitchedMask, mask)

        # Handle overlap based on whether adjustment is necessary
        if (adjust == True) and (i > 0):
            # Adjust image to minimize offset
            img = adjust_image(primaryImg=stitchedImg, secondaryImg=img,
                primaryMask=stitchedMask, secondaryMask=mask,
                expectation=expectation,
                verbose=verbose)

        # Add to stitched image ...
        stitchedImg += img

        # ... and divide by two in the region of overlap
        stitchedImg[overlap == 1] /= 2

        # Update stitched mask
        stitchedMask += mask
        stitchedMask[overlap == 1] /= 2

    return stitchedImg, stitchedMask


def adjust_image(primaryImg, secondaryImg, primaryMask, secondaryMask, expectation='mean', verbose=False):
    '''
    Adjust a given image to seamlessly fit with the reference image.
    Use two categories of methods:
        Differencing finds the difference with the area of overlap
        Bridging fits a 1, 2, or 3-order surface for image that don't overlap
    '''
    # Determine method to use
    if expectation.lower() in ['mean', 'median']:
        # Use differencing method
        adjustedImg = diff_images(primaryImg, secondaryImg,
            primaryMask, secondaryMask,
            expectation,
            verbose)

    elif expectation[:6].lower() == 'bridge':
        # Use bridging method
        adjustedImg = bridge_images(primaryImg, secondaryImg,
            primaryMask, secondaryMask,
            expectation,
            verbose)
    else:
        print('Expectation operator {:s} not recognized'.format(expectation))
        exit()

    return adjustedImg


def diff_images(primaryImg, secondaryImg, primaryMask, secondaryMask, expectation='mean', verbose=False):
    '''
    Adjust the secondary image to the primary image by differencing the two where they overlap.
    '''
    # Explain what is happening
    if verbose == True: print('Adjusting image based on {:s} difference of overlap...'.format(expectation))

    # Find area of overlap
    overlap = find_mask_overlap(primaryMask, secondaryMask)

    # Kill if overlap is zero
    assert np.sum(overlap.flatten()) > 0, 'No overlap, try -e bridge'

    # Compute difference
    D = primaryImg[overlap==1] - secondaryImg[overlap==1]

    # Define expectation operator
    if expectation in ['mean']:
        expct = np.mean
    elif expectation in ['median']:
        expct = np.median

    # Expected difference
    Dexpect = expct(D.flatten())

    # Add difference to image
    secondaryImg += Dexpect

    # Reset masked values to zero
    secondaryImg[secondaryMask == 0] = 0

    # Report if requested
    if verbose == True: 
        print('Expected difference: {:f}'.format(Dexpect))

        # Recompute difference
        Dcorr = primaryImg[overlap==1] - secondaryImg[overlap==1]

        print('Corrected difference: {:f} +- {:f}'.format(expct(Dcorr), np.std(Dcorr)))

    return secondaryImg


def bridge_images(primaryImg, secondaryImg, primaryMask, secondaryMask, expectation='bridge', verbose=False):
    '''
    Fit an nth-degree polynomial surface to the primary image, then adjust the secondary image to that surface.
    '''
    # Explain what is happening
    if verbose == True: print('Adjusting image based on polynomial surface...')

    # Find degree of polynomial
    try:
        # Use last digit of expectation name, e.g., bridge2 = 2nd degree polynomial
        degree = int(expectation[-1])
    except:
        # Assume linear ramp
        degree = 1

    # Report degree if requested
    if verbose == True: print('fitting {:d}-degree polynomial'.format(degree))

    # Fit a surface
    surface, B = fit_surface(primaryImg, primaryMask, degree, decimation=1, verbose=verbose)

    # Subtract surface from secondary image
    D = surface - secondaryImg

    # Mask difference and find expected value of offset
    D = np.ma.array(D, mask=(secondaryMask==0))  # mask difference
    D = D.compressed().flatten()  # compress to 1D array
    Dexpect = np.median(D)  # expected value

    # Remove expected 
    secondaryImg += Dexpect

    # Reset masked values to zero
    secondaryImg[secondaryMask == 0] = 0

    # Residuals of secondary image
    RSS = np.sqrt(np.sum(D**2))
    exptRes = RSS/len(D)

    # Report if requested
    if verbose == True:
        print('Expected difference: {:f}'.format(Dexpect))
        print('Expected residual: {:f}'.format(exptRes))

    return secondaryImg



### PLOTTING ---
def plot_inputs(datasets, masks):
    '''
    Plot input data sets, pairwise. Compute the difference between pairs.
    '''
    # List of data sets
    dsNames = list(datasets.keys())
    nDatasets = len(dsNames)

    # Get map extent
    extent = DS_to_extent(datasets[dsNames[0]])

    # Get primary image and mask
    dsName = dsNames[0]
    primaryImg = datasets[dsName].GetRasterBand(1).ReadAsArray()
    primaryMask = masks[dsName]

    # Loop through each input data set
    for i in range(1, nDatasets):
        # Spawn new figure
        fig, [axPrim, axComp, axDiff] = plt.subplots(ncols=3)
        cbarOrient = 'horizontal'

        # Dataset name
        dsName = dsNames[i]

        # Retrieve images and mask
        compareImg = datasets[dsName].GetRasterBand(1).ReadAsArray()
        compareMask = masks[dsName]

        # Plot primary image
        plot_raster(primaryImg, mask=primaryMask, extent=extent,
            cmap='jet', cbarOrient=cbarOrient, minPct=1, maxPct=99,
            fig=fig, ax=axPrim)

        # Plot comparison image
        plot_raster(compareImg, mask=compareMask, extent=extent,
            cmap='jet', cbarOrient=cbarOrient, minPct=1, maxPct=99,
            fig=fig, ax=axComp)

        # Compute difference between figures
        diffImg = primaryImg - compareImg  # compute strict difference
        diffMask = find_mask_overlap(primaryMask, compareMask)

        plot_raster(diffImg, mask=diffMask, extent=extent,
            cmap='jet', cbarOrient=cbarOrient,
            fig=fig, ax=axDiff)

        # Format figure
        axPrim.set_title('Primary image')
        axComp.set_title('Compare image')
        axDiff.set_title('Difference')


def plot_outputs(stitchedImg, stitchedMask, extent):
    '''
    Plot output data sets.
    '''
    # Spawn figure
    fig, [axImg, axMsk] = plt.subplots(ncols=2)
    cbarOrient = 'vertical'

    # Find clip values
    vmin, vmax = image_percentiles(stitchedImg)

    # Plot stitched image
    plot_raster(stitchedImg, mask=stitchedMask, extent=extent,
        cmap='jet', cbarOrient=cbarOrient, minPct=1, maxPct=99,
        fig=fig, ax=axImg)

    # Plot mask
    plot_raster(stitchedMask, extent=extent,
        cmap='magma', cbarOrient=cbarOrient, vmin=0, vmax=1,
        fig=fig, ax=axMsk)

    # Format plot
    axImg.set_title('Stitched image')
    axMsk.set_title('Stitched mask')



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Load data sets
    if inps.verbose == True: print('*'*32)
    datasets = load_gdal_datasets(inps.imgFiles, verbose=inps.verbose)


    ## Resample to same bounds
    if inps.verbose == True: print('*'*32)
    datasets = match_rasters(datasets, verbose=inps.verbose)


    ## Masking
    if inps.verbose == True: print('*'*32)
    masks = mask_datasets(datasets, inps.maskArgs, verbose=inps.verbose)


    ## Plot formatted inputs if requested
    if inps.plotInputs == True: plot_inputs(datasets, masks)


    ## Stitch data sets
    if inps.verbose == True: print('*'*32)
    stitchedImg, stitchedMask = stitch_datasets(datasets, masks,
        adjust=inps.adjust, expectation=inps.expectation,
        verbose=inps.verbose)


    ## Plot outputs if requested
    if inps.plot == True: plot_outputs(stitchedImg, stitchedMask, DS_to_extent(pick_dataset(datasets)))


    ## Save dataset to file
    # Save as GDAL dataset
    confirm_outdir(inps.outName)
    outName = confirm_outname_ext(inps.outName, ['tif', 'tiff'])
    save_gdal_dataset(outName, stitchedImg, mask=stitchedMask, exDS=pick_dataset(datasets), verbose=inps.verbose)


    plt.show()