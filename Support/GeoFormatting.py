'''
SHORT DESCRIPTION
Geographic coordinate formatting, especially for GDAL data set compatibility.

FUTUTRE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import sys
import numpy as np
from osgeo import gdal
import h5py


### TRANSFORMATIONS ---
class parse_transform:
    def __init__(self, tnsf, M, N, verbose=False):
        '''
        Parse the GDAL transform into full geospatial information.
        '''
        # Parse transform
        self.left, self.dx, self.xshear, self.top, self.yshear, self.dy = tnsf

        # Compute bounds
        self.right = self.left + self.dx*N
        self.bottom = self.top + self.dy*M

        # Record map size
        self.M = M
        self.N = N

        # Report if requested
        if verbose == True:
            print('Spatial extent:')
            print('x-range: {:.3f} - {:.3f}'.format(self.left, self.right))
            print('y-range: {:.3f} - {:.3f}'.format(self.bottom, self.top))
            print('resolution: dx {:.5f}, dy {:.5f}'.format(self.dx, self.dy))


def transform_to_extent(tnsf, M, N, verbose=False):
    '''
    Convert the GDAL transform information to extent for imshow plotting.
    '''
    # Parse spatial info
    spatInfo = parse_transform(tnsf, M, N)

    # Extent
    extent = (spatInfo.left, spatInfo.right, spatInfo.bottom, spatInfo.top)

    # Report if requested
    if verbose == True: print('extent: {:.3f} {:.3f} {:.3f} {:.3f}'.format(*extent))

    return extent


def transform_to_bounds(tnsf, M, N, verbose=False):
    '''
    Convert the GDAL transform information to bounds to gdalwarp resampling.
    '''
    # Parse spatial info
    spatInfo = parse_transform(tnsf, M, N)

    # Bounds
    bounds = (spatInfo.left, spatInfo.bottom, spatInfo.right, spatInfo.top)

    # Report if requested
    if verbose == True: print('bounds: {:.3f} {:.3f} {:.3f} {:.3f}'.format(*bounds))

    return bounds


def xy_to_lola(tnsf, px, py, verbose=False):
    '''
    Convert image coordinates to geo coordinates using transform approach.
     l = T * p + l0
     lon  =  |  dx   yshear| * x + lon0
     lat     |xshear   dy  |   y   lat0
    '''
    # Convert transform to matrix, vector
    T, L0 = vectorize_gdal_transform(tnsf)

    # Formulate position as vector
    p = np.array([[px, py]]).T

    # Transform
    lon, lat = T.dot(p) + L0

    # Flatten results
    lon = lon.flatten()[0]
    lat = lat.flatten()[0]

    # Report if requested
    if verbose == True:
        print('Pixels x/y ({:f}, {:f}) -> Lon/Lat ({:f}, {:f})'.format(px, py, lon, lat))

    return lon, lat


def lola_to_xy(tnsf, lon, lat, verbose=False):
    '''
    Convert geo coordinates to image x,y using transform approach.
     p = Tinv.(l - l0)
     x  =  (|  dx   yshear|)-1 * (lon - lon0)
     y     (|xshear   dy  |)     (lat - lat0)
    '''
    # Convert transform to matrix, vector
    T, L0 = vectorize_gdal_transform(tnsf)

    # Formulate geo coordinates as vector
    L = np.array([[lon, lat]]).T

    # Transform
    Tinv = np.linalg.inv(T)  # inverse transform
    px, py = Tinv.dot((L-L0))

    # Flatten results
    px = int(px.flatten()[0])
    py = int(py.flatten()[0])

    # Report if requested
    if verbose == True:
        print('Lat/Lon ({:f}, {:f}) -> pixels y/x ({:d}, {:d})'.format(lat, lon, py, px))

    return px, py


def vectorize_gdal_transform(tnsf, verbose=False):
    '''
    Format transform matrix and map origin from GDAL geo transform.
     |  dx   yshear|
     |xshear   dy  |
    '''
    # Parse transform data
    left, xstep, xskew, top, yskew, ystep  = tnsf

    # Formulate transform matrix
    T = np.array([[xstep, yskew],
                  [xskew, ystep]])

    # Origin points as vector
    L0 = np.array([[left, top]]).T

    # Show if requested
    if verbose == True:
        print('T\t|xstep {:.4f}, yskew {:.4f}|\n \t|xskew {:.4f}, ystep {:.4f}|'.format(xstep, yskew, xskew, ystep))
        print('Origin: |left: {:.4f}|\n        | top: {:.4}|'.format(left, top))

    return T, L0



### GDAL FORMATTING---
def get_raster_size(DS):
    '''
    Retrieve raster size from data set.
    '''
    M = DS.RasterYSize
    N = DS.RasterXSize

    return M, N


def DS_to_extent(DS, verbose=False):
    '''
    Extract the geographic extent in pyplot imshow format from a GDAL data set.
    '''
    # Raster size
    M, N = get_raster_size(DS)

    # Get transform
    tnsf = DS.GetGeoTransform()

    # Convert to extent
    extent = transform_to_extent(tnsf, M, N, verbose=False)

    return extent


def DS_to_bounds(DS, verbose=False):
    '''
    Extract the geographic bounds in GDAL warp format form a GDAL data set.
    '''
    # Raster size
    M, N = get_raster_size(DS)

    # Get transform
    tnsf = DS.GetGeoTransform()

    # Convert to bounds
    bounds = transform_to_bounds(tnsf, M, N, verbose=False)

    return bounds


def determine_common_bounds(datasets, cropping='union', resolution='fine', verbose=False):
    '''
    Determine the common bounds based on a list of provided GDAL data sets.
    INPUTS
        datasets is a list of GDAL data sets
        cropping determines how to treat the overlap of datasets ([union], intersection)
        resolution gives the finest or coarsest resolution, based on the inputs ([fine], coarse)
    OUTPUTS
        bounds are the map extent in GDAL format(xmin, ymin, xmax, ymax)
        xRes, yRes are the pixel sizes
    '''
    # Setup
    bounds = []  # empty bounds
    xRes = []  # empty pixel x-size
    yRes = []  # empty pixel y-size

    # Loop through data sets
    for dataset in datasets:
        # Retrieve geotransform
        tnsf = dataset.GetGeoTransform()

        # Retrieve raster size
        M, N = dataset.RasterYSize, dataset.RasterXSize

        # Retrieve bounds
        bounds.append(transform_to_bounds(tnsf, M, N, verbose=verbose))

        # Retrieve pixel resolution
        _, dx, _, _, _, dy = tnsf
        xRes.append(dx)
        yRes.append(np.abs(dy))

    # Format bounds
    left = [bound[0] for bound in bounds]
    bottom = [bound[1] for bound in bounds]
    right = [bound[2] for bound in bounds]
    top = [bound[3] for bound in bounds]

    # Determine bounds
    if cropping == 'union':
        xmin = min(left)
        ymin = min(bottom)
        xmax = max(right)
        ymax = max(top)
    elif cropping == 'intersection':
        xmin = max(left)
        ymin = max(bottom)
        xmax = min(right)
        ymax = min(top)

    bounds = (xmin, ymin, xmax, ymax)

    # Determine resolution
    if resolution == 'fine':
        xRes = min(xRes)
        yRes = min(yRes)
    elif resolution == 'coarse':
        xRes = max(xRes)
        yRes = max(yRes)

    # Report if requested
    if verbose == True:
        print('Global bounds: {:.3f} {:.3f} {:.3f} {:.3f}'.format(*bounds))
        print('Resolution: {:.5f} x {:.5f}'.format(xRes, yRes))

    return bounds, xRes, yRes


def grid_from_DS(DS, verbose=False):
    '''
    Create an X, Y grid given a GDAL data set.
    '''
    if verbose == True: print('Creating grid from data set')

    # Get map size
    M, N = get_raster_size(DS)

    # Retrieve geographic transform
    tnsf = DS.GetGeoTransform()

    # Parse geographic information
    spatInfo = parse_transform(tnsf, M, N)

    # Create grid
    X, Y = grid_from_spatial_info(spatInfo)

    return X, Y



### MINTPY FORMATTING ---
def get_mintpy_reference_point(DSname, verbose=False):
    '''
    Retrieve the reference point in lon/lat from a MintPy data set.
    '''
    # Load attributes
    with h5py.File(DSname, 'r') as DS:
        refLon = float(DS.attrs['REF_LON'][:])
        refLat = float(DS.attrs['REF_LAT'][:])

    # Report if requested
    if verbose == True:
        print('Reference point: Lon {:.6f}, Lat {:.6f}'.format(refLon, refLat))

    return refLon, refLat


def get_mintpy_transform(DSname, verbose=False):
    '''
    Retrieve and parse the geographic metadata of a MintPy data set.

    INPUTS
        DSname is the name of a MintPy data set
    '''
    # Load attributes from MintPy data set
    with h5py.File(DSname, 'r') as DS:
        attrs = DS.attrs

        # Parse attributes
        tnsf = parse_mintpy_geo_attributes(attrs, verbose=verbose)

    return tnsf


def parse_mintpy_geo_attributes(attrs, verbose=False):
    '''
    Parse the geographic attributes of a MintPy data set.
    '''
    # Parse geographic attributes
    LonRefs = [attrs['LON_REF1'], attrs['LON_REF2'], attrs['LON_REF3'], attrs['LON_REF4']]
    LatRefs = [attrs['LAT_REF1'], attrs['LAT_REF2'], attrs['LAT_REF3'], attrs['LAT_REF4']]
    xstep = attrs['X_STEP']
    ystep = attrs['Y_STEP']

    # Convert to floats
    LonRefs = [float(lon) for lon in LonRefs]
    LatRefs = [float(lat) for lat in LatRefs]
    xstep = float(xstep)
    ystep = float(ystep)

    left = np.min(LonRefs)
    top = np.max(LatRefs)

    # Find transform values
    tnsf = [left, xstep, 0, top, 0, ystep]

    # Report if requested
    if verbose == True:
        print('MintPy geo transform: {:.8f} {:.6f} {:.4f} {:.8f} {:.4f} {:.6f}'.format(*tnsf))

    return tnsf



### GRIDDING ---
def grid_from_spatial_info(spatInfo, verbose=False):
    '''
    Create an X, Y grid from a spatInfo object.
    '''
    if verbose == True: print('Creating spatial grid')

    # Create grid from spatial information object
    x = np.linspace(spatInfo.left, spatInfo.right, spatInfo.N)
    y = np.linspace(spatInfo.top, spatInfo.bottom, spatInfo.M)

    X, Y = np.meshgrid(x, y)

    return X, Y


def grid_from_transform(tnsf, M, N, verbose=False):
    '''
    Create and X, Y grid from a geographic transform and map size.
    '''
    if verbose == True: print('Creating spatial grid')

    # Formulate spatial information
    spatInfo = parse_transform(tnsf, M, N)

    # Create X, Y grid from spatial information object
    X, Y = grid_from_spatial_info(spatInfo)

    return X, Y


def rotate_grid(X, Y, angle):
    '''
    '''
    pass