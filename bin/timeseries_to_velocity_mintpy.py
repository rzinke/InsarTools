#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Convert a MintPy displacement timeseries to velocities.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py
from IOsupport import confirm_outdir, confirm_outname_ext, append_fname, load_mintpy_timeseries, save_gdal_dataset
from GeoFormatting import get_mintpy_transform, transform_to_extent, lola_to_xy
from Masking import create_mask
from Viewing import plot_raster
from Fitting import dates_to_datetimes, time_since_reference, velocity_from_timeseries


### PARSER ---
Description = '''Convert a MintPy displacement timeseries to velocities.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='tsName', type=str,
        help='MintPy timeseries.h5 file.')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')
    InputArgs.add_argument('-r','--reference-point', dest='refPoint', type=float, nargs=2, default=None,
        help='Reference point ([None]; lon lat).')
    InputArgs.add_argument('-f','--fit-type', dest='fitType', type=str, default='linear',
        help='Fit type used to compute velocities.')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true', 
        help='Plot.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### REFERENCING ---
def adjust_to_refpoint(disps, tnsf, refLon, refLat, verbose=False):
    '''
    Calculate the reference point value in pixel coordinates and subtract it
     from each image in the timeseries.
    '''
    if verbose == True: print('Removing reference point: Lon {:f} Lat {:f}'.format(refLon, refLat))

    # Convert geographic coordinates to pixels
    px, py = lola_to_xy(tnsf, refLon, refLat, verbose=verbose)

    # Remove reference pixel value from each image
    for i in range(disps.shape[0]):
        disps[i,:,:] = disps[i,:,:]-disps[i,py,px]

    return disps



### SAVING ---
def save_velocities(outName, velocities, mask, tnsf, verbose):
    '''
    Save each component of the velocity solution.
    '''
    if verbose == True: print('Saving results')

    # Format output name
    outName = confirm_outname_ext(outName, ext=['tif'])
    confirm_outdir(outName)

    # Retrieve components list
    componentsList = velocities.componentsList

    # Loop through map components
    for component in componentsList:
        # Format specific output name
        fname = append_fname(outName, '_{:s}'.format(component.replace('Map', '')))

        # Save to GDAL data set
        save_gdal_dataset(fname, velocities.__dict__[component],
            mask=mask, proj='WGS84', tnsf=tnsf,
            verbose=verbose)



### PLOTTING ---
def plot_velocities(velocities, mask, extent, verbose=False):
    '''
    Plot each component of the velocity map.
    '''
    if verbose == True: print('Plotting results')

    # Retrieve components list
    componentsList = velocities.componentsList

    # Loop through map components
    for component in componentsList:
        if verbose == True:
            print(component)
            fig, ax = plot_raster(velocities.__dict__[component], mask=mask, extent=extent,
                cbarOrient='auto', minPct=1, maxPct=99)
            ax.set_title(component.replace('Map', ''))



### MAIN ---
if __name__ == '__main__':
    ## Setup
    # Gather arguments
    inps = cmdParser()


    ## Load data
    # Open and close HDF5 data set
    dates, disps = load_mintpy_timeseries(inps.tsName, verbose=inps.verbose)

    # Convert dates to datetimes
    datetimes = dates_to_datetimes(dates, verbose=inps.verbose)

    # Calculate time since beginning of series
    times = np.array(time_since_reference(datetimes, verbose=inps.verbose))

    # Create mask
    mask = create_mask(disps[-1,:,:], inps.maskArgs, verbose=inps.verbose)

    # Data set sizes
    nEpochs, M, N = disps.shape

    # Geographic transform
    tnsf = get_mintpy_transform(inps.tsName, verbose=inps.verbose)

    # Geographic extent
    extent = transform_to_extent(tnsf, M, N, verbose=inps.verbose)


    ## Referencing
    # Remove reference point if specified
    if inps.refPoint is not None:
        # Adjust to reference point
        disps = adjust_to_refpoint(disps, tnsf, inps.refPoint[0], inps.refPoint[1], verbose=inps.verbose)


    ## Timeseries to velocity
    velocities = velocity_from_timeseries(disps, times, mask=mask, fitType=inps.fitType, verbose=inps.verbose)


    ## Save results
    save_velocities(outName=inps.outName, velocities=velocities, mask=mask, tnsf=tnsf, verbose=inps.verbose)


    ## Plot results
    if inps.plot == True:
        plot_velocities(velocities=velocities, mask=mask, extent=extent, verbose=inps.verbose)

        plt.show()