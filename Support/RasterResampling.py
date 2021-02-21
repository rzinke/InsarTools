'''
SHORT DESCRIPTION
Resample rasters as GDAL data sets.

INHERITANCES
GeoFormatting: determine_common_bounds

TESTING STATUS
match_rasters 'intersection' feature needs to be tested.
'''

### IMPORT MODULES ---
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
        datasets is a list or dictionary of resampled GDAL data sets
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