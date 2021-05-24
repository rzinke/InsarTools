'''
SHORT DESCRIPTION
Resample rasters as GDAL data sets.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import numpy as np
from osgeo import gdal
from GeoFormatting import determine_common_bounds, get_raster_size


### RASTER RESAMPLING---
def gdal_resample(dataset, bounds, xRes=None, yRes=None, M=None, N=None, alg='near',
    outName=None, verbose=False):
    '''
    Resample data set using GDAL warp.
    '''
    # Setup outputs
    if outName is None:
        outName = ''
        fmt = 'MEM'
    else:
        fmt = 'GTiff'

    # Resample
    dataset = gdal.Warp(outName, dataset,
        options=gdal.WarpOptions(format=fmt, outputBounds=bounds,
            xRes=xRes, yRes=yRes, width=N, height=M,
            resampleAlg=alg))

    return dataset


def match_rasters(datasets, cropping='union', resolution='fine', verbose=False):
    '''
    Resample rasters to a common grid.
    INPUTS
        datasets is a list or dictionary of georeferenced GDAL data sets
        cropping determines how to treat the overlap of datasets ([union], intersection)
        resolution gives the finest or coarsest resolution, based on the inputs ([fine], coarse)
    OUTPUTS
        datasets is a dictionary or list of resampled GDAL data sets
    '''
    if verbose == True: print('Resampling rasters to common bounds')

    # Check that datasets is a list
    dsNames = None  # assume no names
    if type(datasets) == dict:
        dsNames = list(datasets.keys())
        datasets = list(datasets.values())

    # Size of original raster
    M, N = get_raster_size(datasets[0])

    # Determine common bounds and resolution
    bounds, xRes, yRes = determine_common_bounds(datasets, cropping, resolution, verbose)

    # Loop through to resample data sets
    resDatasets = []  # empty list to store resampled data sets
    for dataset in datasets:
        # Resample data set
        resDatasets.append(gdal_resample(dataset, bounds, xRes=xRes, yRes=yRes))

    # Reformat as dictionary if input type was dictionary
    if dsNames is not None:
        resDatasets = dict(zip(dsNames, resDatasets))

    return resDatasets



### POINTS IN RASTER ---
def sample_points_from_raster(img, x, y, mask=None, searchR=0, verbose=False):
    '''
    Sample points with the given coordinates from the provided raster.
    All units are given in pixel coordinates.
    '''
    if verbose == True: print('Sampling points from raster')

    # Parameters
    nPts = len(x)  # number of points to sample
    searchR = int(searchR)  # ensure search radius is an integer
    w = np.arange(-searchR, searchR+1)  # width of search kernel (pixels)

    if verbose == True:
        print('\t{:d} sample points'.format(nPts))
        print('\tusing search radius of {:d} pixels'.format(searchR))

    # Check mask
    if mask is None: mask = np.ones(img.shape)

    # Loop through sample points
    sX =  []  # valid sample x coordinates
    sY =  []  # valid sample y coordinates
    sZ =  []  # valid sample values
    indices = []  # indices of valid samples

    for i in range(nPts):
        # Gather image data within search radius
        data = img[y[i]+w, x[i]+w]

        # Exclude masked pixels
        localMask = mask[y[i]+w, x[i]+w]
        data = data[localMask == 1]  # use only unmasked data

        # Record valid sample values
        if data.size > 0:
            # Expected value
            expcVal = np.median(data)

            # Append to list
            sX.append(x[i])
            sY.append(y[i])
            sZ.append(expcVal)
            indices.append(i)

    # Convert data to numpy arrays
    sX = np.array(sX, dtype=int)
    sY = np.array(sY, dtype=int)
    sZ = np.array(sZ)

    # Convert indices to index array
    ndx = [False]*nPts
    for i in indices:
        ndx[i] = True

    # Update stats
    nPts = len(sZ)

    # Report if requested
    if verbose == True: print('\tsampled {:d} valid points'.format(nPts))

    return sX, sY, sZ, ndx