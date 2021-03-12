'''
SHORT DESCRIPTION
Data set loading and saving, especially for GDAL compatibility.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import os
import numpy as np
from osgeo import gdal
import h5py


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


def parse_extension(fname, verbose=False):
    '''
    Retrieve the filename exension.
    '''
    # Get extension
    ext = fname.split('.')[-1]

    # Report if requested
    if verbose == True: print('Filename extension: {:s}'.format(ext))

    return ext


def confirm_overwrite(fname):
    '''
    Check if file already exists and return True to confirm overwrite.
    '''
    # Check if file already exists
    if os.path.exists(fname):
        # Check file existence
        response = input('File already exists. Overwrite? (y/n) ')

        # Record user response
        if response.lower() in ['y', 'yes']:
            write_file = True  # overwrite
        else:
            write_file = False  # don't overwrite
    else:
        # File does not already exist
        write_file = True  # write new file

    return write_file



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
    check_exists(dsPath)

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
        ourName - data set output name
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



### HDF5 LOADING ---
def load_hdf5_dataset(dsName, verbose=False):
    '''
    Check that the specified HDF5 data set exists and load it.
    '''
    # Check that HDF5 file exists
    check_exists(dsName)

    # Load HDF5 data set
    DS = h5py.File(dsName, 'r')

    # Report if requested
    if verbose == True: print('Loaded HDF5 file: {:s}'.format(dsName))

    return DS


def load_mintpy_velocity(velName, verbose=False):
    '''
    Load and parse a MintPy velocity.h5 file.
    '''
    # Load file after checking it exists
    with load_hdf5_dataset(velName, verbose=verbose) as DS:
        # Parse file contents
        velocity = DS['velocity'][:]  # velocity

    # Report if requested
    if verbose == True:
        print('Parsed MintPy velocity: {:d} x {:d} map'.format(*velocity.shape))

    return velocity


def load_mintpy_timeseries(tsName, verbose=False):
    '''
    Load and parse a MintPy timeseries.h5 file.
    '''
    # Load file after checking it exists
    with load_hdf5_dataset(tsName, verbose=verbose) as DS:
        # Parse file contents
        dates, disps = parse_mintpy_timeseries(DS, verbose=verbose)

    return dates, disps


def parse_mintpy_timeseries(DS, verbose=False):
    '''
    Convenience function for parsing a MintPy timeseries data set.
    '''
    # Parse data set
    dates = [date.astype(str) for date in DS['date']]  # dates
    disps = DS['timeseries'][:]  # timeseries

    # Report if requested
    if verbose == True:
        print('Parsed MintPy displacements and dates')
        print('{:d} dates'.format(len(dates)))
        print('{:d} displacements, in {:d} x {:d} maps'.format(*disps.shape))

    return dates, disps


def load_mintpy_geometry(geomName, verbose=False):
    '''
    Load MintPy geometry data set.
    '''
    # Load file after checking it exists
    with load_hdf5_dataset(geomName, verbose=verbose) as DS:
        # Parse file contents
        inc = DS['incidenceAngle'][:]  # incidence
        az = DS['azimuthAngle'][:]  # azimuth

    # Report if requested
    if verbose == True:
        print('Loaded MintPy geometry file')
        print('Parsed incidence and azimuth')
        print('Sizes: {:d} x {:d}'.format(*inc.shape))



### POLYLINE LOADING ---
def load_polyline(polylineName, verbose=False):
    '''
    Load a polyline from a text file in x, y format.
    '''
    if verbose == True: print('Loading polyline file...')

    # Check file exists
    check_exists(polylineName)

    # Load file contents
    data = np.loadtxt(polylineName)

    # Parse contents
    x = data[:,0]
    y = data[:,1]

    nVertices = data.shape[0]

    # Report if requested
    if verbose == True:
        print('{:d} vertices detected'.format(nVertices))

    return x, y



### MISCELLANEOUS ---
def check_exists(fname):
    '''
    Check that a file exists.
    '''
    assert os.path.exists(fname), 'ERROR: {:s} does not exist!!!'.format(fname)


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