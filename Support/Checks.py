'''
SHORT DESCRIPTION
Sanity checks.

FUTUTRE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import sys
import matplotlib.pyplot as plt
from GeoFormatting import get_raster_size, parse_transform, lola_to_xy
from Viewing import *
import time


### DATASET FORMATTING ---
def check_dataset_sizes(datasets):
    '''
    Check that all data sets and masks have the same dimensions in pixels.
    '''
    # Setup
    dsNames = list(datasets.keys())
    nDatasets = len(dsNames)

    # Initial data set size
    M0, N0 = get_raster_size(datasets[dsNames[0]])

    # Loop through data sets
    for i in range(1, nDatasets):
        # Compare data set sizes
        dsName = dsNames[i]
        M, N = get_raster_size(datasets[dsName])

        assert (M0, N0) == (M, N), \
            'Dimensions of data set {:s} ({:d} x {:d}) not consistent with those of reference data set ({:d} x {:d}).'.\
            format(dsName, M, N, M0, N0)

    return M0, N0


def check_loc_in_dataset(DS, mask, lon, lat, verbose=False):
    '''
    Check that the lon/lat are within the bounds of a dataset, and that the value is not masked.
    '''
    if verbose == True: print('Checking point location')

    # Retrieve map information
    tnsf = DS.GetGeoTransform()  # geo transform
    M, N = get_raster_size(DS)  # raster dimensions
    geos = parse_transform(tnsf, M, N)  # transform class

    # Check that given coordinates are within the confines of the data set
    if (lon < geos.left) or (lon > geos.right) or (lat < geos.bottom) or (lat > geos.top):
        print('Error: Lon/Lat ({:f}E, {:f}N) outside of map extent ({left:f}-{right:f}E, {bottom:f}-{top:f}N)'.\
            format(lon, lat, **geos.__dict__), file=sys.stderr)
        exit()

    # Convert lon/lat to pixel location
    px, py = lola_to_xy(tnsf, lon, lat, verbose=verbose)

    # Check that the pixels are do not coincide with a masked value
    if mask[py, px] == 0:
        # Show where the pixel is located, within the mask
        img = DS.GetRasterBand(1).ReadAsArray()
        fig, axes = plot_map_mask(img, mask, extent=transform_to_extent(tnsf, M, N))  # plot data set
        [axis.scatter(lon, lat, facecolor='w', edgecolor='k') for axis in axes]  # plot point location
        fig.suptitle('Error. Pixel location falls within a masked value')

        # Print error message
        print('Error: Pixel (x {:d}, y {:d}) falls within a masked value'.format(px, py),
            file=sys.stderr)

        plt.show()
        exit()

    # Report if requested
    if verbose == True: print('Checks passed.')

    return px, py



### MISCELLANEOUS ---
def timed_run(fcn, label=None, verbose=False):
    '''
    Clock the time it takes to complete a function.
    '''
    # Start clock
    startTime = time.time()

    # Run function
    outputs = fcn

    # End clock
    endTime = time.time()

    # Report if requested
    if verbose == True: print('{:f} to complete'.format(endTime-startTime), label)

    return outputs