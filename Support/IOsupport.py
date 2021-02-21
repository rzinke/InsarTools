'''
SHORT DESCRIPTION
Data set loading and saving, especially for GDAL compatibility.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import os
from osgeo import gdal


### NAME FORMATTING ---
def rel_to_abs_paths(fnames):
    '''
    Convert relative file paths to absolute file paths.
    '''
    # Use list if not already that type
    if type(fnames) != list: fnames = [fnames]

    # Convert relative paths to absolute paths
    fnames = [os.path.abspath(fname) for fname in fnames]

    return fnames


def confirm_outdir(fname):
    '''
    Confirm that the output folder exists.
    '''
    # Check that output name is an absolute file path
    fname = os.path.abspath(fname)

    # Get directory
    dirName = os.path.dirname(fname)

    # Confirm directory exists
    if not os.path.exists(dirName):
        os.mkdir(dirName)


def append_fname(fname, appendix, verbose=False):
    '''
    "Append" a filename by placing the appendix after the main file name and
     before the extension.
    '''
    parts = fname.split('.')  # split extension
    ext = parts[-1]  # extension
    parts = '.'.join(parts[:-1])  # rejoin non-extension parts
    parts += appendix  # add appendix
    appName = '.'.join([parts, ext])  # reform with extension

    if verbose == True: print('Filename {:s} -> {:s}'.format(fname, parts))

    return appName


def confirm_outname_ext(outName, ext=['tif', 'tiff'], verbose=False):
    '''
    Check that the output name uses the specified extension(s). If not, add that extension.
    '''
    # Check that extension is used
    if outName.split('.')[-1] not in ext:
        outName = '.'.join([outName, ext[0]])

    # Report if reequested
    if verbose == True: print('Out name: {:s}'.format(outName))

    return outName



### LOADING GDAL DATASETS ---
def load_gdal_datasets(dsPaths, dsNames=None, verbose=False):
    '''
    Provide a list of dataset names.
    Check that the files actually exist, then load into dictionary.
    '''
    if verbose == True: print('Loading multiple data sets')

    # Create list of data set names if not already provided
    if dsNames is None:
        dsNames = []  # empty list of names
        for i, dsPath in enumerate(dsPaths):
            # Formulate name
            dsName = os.path.basename(dsPath)  # path basename
            dsName = dsName.split('.')[:-1]  # remove extension
            dsName = '.'.join(dsName)  # reformulate name

            # Check that name not already in list
            if dsName in dsNames: dsName += '_{:d}'.format(i)

            dsNames.append(dsName)

    # Setup
    datasets = {}

    # Loop through data set names
    for i, dsPath in enumerate(dsPaths):
        # Data set name
        dsName = dsNames[i]

        # Report if requested
        if verbose == True: print('Loading {:s}: {:s}'.format(dsName, dsPath))

        # Load data set
        datasets[dsName] = load_gdal_dataset(dsPath)

    return datasets


def load_gdal_dataset(dsPath, verbose=False):
    '''
    Load data set using GDAL.
    '''
    if verbose == True: print('Loading data set')

    # Check that data set path exists
    if not os.path.exists(dsPath):
        print('ERROR: {:s} does not exist!!!'.format(dsPath))
        exit()

    # Open GDAL data set
    DS = gdal.Open(dsPath, gdal.GA_ReadOnly)

    # Report if requested
    if verbose == True: print('Loaded: {:s}'.format(dsPath))

    return DS



### SAVING GDAL DATASETS ---
def save_gdal_dataset(outName, imgs, mask=None, exDS=None, proj=None, tnsf=None, fmt='GTiff', verbose=False):
    '''
    Save image to georeferenced data set, given an example data set.

    INPUTS
        imgs - single image or list of images; if a list is provided, each item
         in the list will be written to a separate band in the order in which they occcur
        (mask) is the mask to be applied to all bands. Masked values set to zero.
        (exDS) provides the projection and geotransform unless they are explicitly
         specified
        (proj) is the geographic projection, overrides exDS
        (tnsf) is the geographic transform, overrides exDS
        (fmt) is the output image format
    '''
    # Spatial parameters
    if proj is None:
        assert exDS is not None, 'Geographic projection or example data set must be provided'
        proj = exDS.GetProjection()
    if tnsf is None:
        assert exDS is not None, 'Geographic transform or example data set must be provided'
        tnsf = exDS.GetGeoTransform()

    # Convert single band to list if not already
    if type(imgs) != list: imgs = [imgs]

    # Image parameters
    nBands = len(imgs)

    # Check that all images are the same size
    M, N = imgs[0].shape
    for img in imgs:
        assert img.shape == (M, N), 'Images must be the same size'

    # Mask image(s) if mask is provided
    if mask is not None: 
        for img in imgs:
            img[mask==0] = 0  # apply mask value of zero

    # Establish driver
    driver = gdal.GetDriverByName(fmt)

    # Create data set
    DSout = driver.Create(outName, N, M, nBands, gdal.GDT_Float32)
    DSout.SetProjection(proj)
    DSout.SetGeoTransform(tnsf)

    # Write image bands
    for n in range(nBands):
        DSout.GetRasterBand(n+1).WriteArray(imgs[n])

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