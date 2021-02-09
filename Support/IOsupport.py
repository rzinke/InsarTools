'''
SHORT DESCRIPTION
Data set loading and saving, especially for GDAL compatibility.

INHERITANCES

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import os
from osgeo import gdal


### LOADING GDAL DATASETS ---
def load_gdal_datasets(dsPaths, dsNames=None, verbose=False):
    '''
    Provide a list of dataset names.
    Check that the files actually exist, then load into dictionary.
    '''
    # Setup
    datasets = {}

    # Loop through data set names
    for i, dsPath in enumerate(dsPaths):
        # Format data set name
        if dsNames is None:
            # Create dataset name based on file name
            dsName = os.path.basename(dsPath)  # path basename
            dsName = dsName.split('.')[:-1]  # remove extension
            dsName = '.'.join(dsName)  # reformulate name
        else:
            # Use existing data set name from list
            dsName = dsNames[i]

        # Report if requested
        if verbose == True: print('Loading: {:s}'.format(dsName))

        # Load data set
        datasets[dsName] = load_gdal_dataset(dsPath, verbose=verbose)

    return datasets


def load_gdal_dataset(dsPath, verbose=False):
    '''
    Load data set using GDAL.
    '''
    # Check that data set path exists
    if not os.path.exists(dsPath):
        print('ERROR: {:s} does not exist!!!'.format(dsPath))
        exit()

    # Open GDAL data set
    DS = gdal.Open(dsPath, gdal.GA_ReadOnly)

    return DS



### SAVING GDAL DATASETS ---
def save_gdal_dataset(outName, img, exDS, fmt='GTiff', verbose=False):
    '''
    Save image to georeferenced data set, given an example data set.
    '''
    # Establish driver
    driver = gdal.GetDriverByName(fmt)

    # Set parameters
    M, N = img.shape
    proj = exDS.GetProjection()
    tnsf = exDS.GetGeoTransform()

    # Create data set
    DSout = driver.Create(outName, N, M, 1, gdal.GDT_Float32)
    DSout.SetProjection(proj)
    DSout.SetGeoTransform(tnsf)
    DSout.GetRasterBand(1).WriteArray(img)

    # Flush cache
    DSout.FlushCache

    # Report if requested
    if verbose == True: print('Saved raster to: {:s}'.format(outName))



### MISCELLANEOUS ---
def pick_dataset(datasets, verbose=False):
    '''
    Return the first entry in a dictionary.
    '''
    # Create list of datasets
    dsNames = list(datasets.keys())

    # Return first dataset
    dsName = dsNames[0]

    # Report if requested
    if verbose == True: print('Returning dictionary entry: {:s}'.format(dsName))

    # Return first data set
    return datasets[dsName]


def check_outname_ext(outName, ext=['.tif', '.tiff'], verbose=False):
    '''
    Check that the output name uses the specified extension(s). If not, add that extension.
    '''
    # Check that extension is used
    if outName.split('.')[-1] not in ext:
        outName += ext[0]

    # Report if reequested
    if verbose == True: print('Out name: {:s}'.format(outName))

    return outName