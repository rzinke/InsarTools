#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Convert timeseries displacements to velocities.

FUTURE IMPROVEMENTS
    * Accept list of fnames
    * Accept list of dates in command line
    * Spatial reference point

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from pandas.plotting import register_matplotlib_converters
import multiprocessing as mp
from IOsupport import confirm_outdir, confirm_outname_ext, append_fname, load_gdal_dataset, save_gdal_dataset
from GeoFormatting import transform_to_extent, lola_to_xy
from Masking import create_mask
from Viewing import plot_raster
from Fitting import dates_to_datetimes, time_since_reference, velocity_from_timeseries


### PARSER ---
Description = '''Convert timeseries displacements to velocities.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='tsName', type=str,
        help='Timeseries file(s).')
    InputArgs.add_argument('-d','--date-list', dest='dateList', type=str, required=True,
        help='List of dates in format YYYYMMDD (list or file with one date per line).')
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
        help='Plot mode.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### LOADING AND FORMATTING ---
def format_displacements(DS, verbose=False):
    '''
    Load images and parse the spatial data.
    '''
    # Number of epochs
    nEpochs = DS.RasterCount

    # Format displacement maps into array
    disps = []
    for i in range(nEpochs):
        disps.append(DS.GetRasterBand(i+1).ReadAsArray())

        # Replace NaNs with zeros
        disps[i][np.isnan(disps[i])==1] = 0

    disps = np.array(disps)

    # Report if requested
    if verbose == True:
        print('Number of epochs detected: {:d}'.format(nEpochs))
        print('Image size: {:d} x {:d}'.format(DS.RasterYSize, DS.RasterXSize))

    return disps


def load_dates(dateList, verbose=False):
    '''
    Load dates from list file.
    Compute the cumulative time since the start date.
    '''
    if verbose == True: print('Loading dates')

    # Load dates from date list
    with open(dateList, 'r') as dateFile:
        dates = dateFile.readlines()

    # Format dates
    dates = [date.strip(',').strip('\n') for date in dates]

    # Convert dates to datetimes
    datetimes = dates_to_datetimes(dates, verbose=verbose)

    # Calculate time since beginning of series
    times = np.array(time_since_reference(datetimes, verbose=verbose))

    return dates, times



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
def save_velocities(outName, velocities, mask, exDS, verbose):
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
        save_gdal_dataset(fname, velocities.__dict__[component].copy(),
            mask=mask, exDS=exDS,
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
            fig, ax = plot_raster(velocities.__dict__[component], mask=mask, extent=extent,
                cbarOrient='auto', minPct=1, maxPct=99)
            ax.set_title(component.replace('Map', ''))



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather inputs
    inps = cmdParser()


    ## Load and format data
    # Load GDAL data set
    DS = load_gdal_dataset(inps.tsName, verbose=inps.verbose)

    # Format epoch displacements
    disps = format_displacements(DS, verbose=inps.verbose)

    # Load and format dates
    dates, times = load_dates(inps.dateList, verbose=inps.verbose)

    # Create mask
    mask = create_mask(disps[-1,:,:], inps.maskArgs, verbose=inps.verbose)

    # Geographic transform
    tnsf = DS.GetGeoTransform()

    # Geographic extent
    extent = transform_to_extent(tnsf, DS.RasterYSize, DS.RasterXSize)


    ## Referencing
    # Remove reference point if specified
    if inps.refPoint is not None:
        # Adjust to reference point
        disps = adjust_to_refpoint(disps, tnsf, inps.refPoint[0], inps.refPoint[1], verbose=inps.verbose)


    ## Timeseries to velocity
    velocities = velocity_from_timeseries(disps, times, mask=mask, fitType=inps.fitType, verbose=inps.verbose)
    print(velocities.velocityMap)


    ## Save results
    save_velocities(outName=inps.outName, velocities=velocities, mask=mask, exDS=DS, verbose=inps.verbose)


    ## Plot results
    if inps.plot == True:
        plot_velocities(velocities=velocities, mask=mask, extent=extent, verbose=inps.verbose)

        plt.show()