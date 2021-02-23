'''
SHORT DESCRIPTION
Conversion factors and routines for different radar geometry conventions.

FUTUTRE IMPROVEMENTS
    * Weighting for inversion scheme
    * Parallelize inversion scheme

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import numpy as np
from IOsupport import load_gdal_dataset
import timeit


### ARIA GEOMETRIES ---
class aria_sensitivity_factors:
    def __init__(self, inc, az, verbose=False):
        '''
        Derive scaling factors for target-to-sensor look vectors and line of sight
         projections based on the incidence and azimuth angles in the ARIA
         convention.
        This follows the ARIA convention:
         https://github.com/aria-tools/ARIA-tools-docs/blob/master/JupyterDocs/ariaExtract/ariaExtract_tutorial.ipynb

        INPUTS
            alpha is the azimuth angle, measured counterclockwise between the east
             and the look direction from the ground target to the sensor
            theta is the incidence angle with respect to vertical
        '''
        # Convert degrees to radians
        inc = np.deg2rad(inc)
        az = np.deg2rad(az)

        # Incidence-based components
        self.vFactor = np.cos(inc)  # vertical
        self.hFactor = np.sin(inc)  # horizontal

        # Horizontal components
        self.eFactor = np.cos(az)  # east
        self.nFactor = np.sin(az)  # north


def aria_geom_to_vector(inc, az, verbose=False):
    '''
    Convert the incidence and azimuth angles to vectors based on the ARIA
     convention and scaling factors given by the aria_sensitivity_factors class.

    INPUTS
        inc and az are from the respective ARIA files, and are given in degrees
    '''
    # Determine scaling factors
    S = aria_sensitivity_factors(inc, az)

    # Formulate pointing vector
    Px = S.hFactor*S.eFactor  # east component
    Py = S.hFactor*S.nFactor  # north component
    Pz = S.vFactor  # vertical component

    # Report if requested
    if verbose == True:
        print('Median pointing vector:')
        print('\tPx: {:.6f}'.format(np.median(Px)))
        print('\tPy: {:.6f}'.format(np.median(Py)))
        print('\tPz: {:.6f}'.format(np.median(Pz)))

    return Px, Py, Pz


def aria_to_los(ux, uy, uz, inc, az, verbose=False):
    '''
    Project 3D displacements into sensor line of sight based on the ARIA 
     convention and scaling factors given by the aria_sensitivity_factors class.
    '''
    # Determine scaling factors
    S = aria_sensitivity_factors(inc, az)

    # Compute line of sight displacement
    LOS = S.hFactor*(ux*S.eFactor + uy*S.nFactor) + uz*S.vFactor

    # Report if requested
    if verbose == True:
        print('Median LOS value: {:.6f}'.format(np.median(LOS)))

    return LOS



### ISCE GEOMETRIES ---
class isce_sensitivity_factors:
    def __init__(self, inc, az, verbose=False):
        '''
        Derive scaling factors for target-to-sensor look vectors and line of sight
         projections based on the incidence and azimuth angles in the ISCE
         convention.

        From the los.rdr.geo.xml file:
        Two channel Line-Of-Sight geometry image (all angles in degrees).
        Represents vector drawn from target to platform.
        Channel 1: Incidence angle measured from vertical at target (always +ve).
        Channel 2: Azimuth angle measured from North in Anti-clockwise direction.

        INPUTS
            * All angles given in degrees
            inc is the incidence angle between the vertical and the LOS
            az is the azimuth angle, measured anti-clockwise from north
        '''
        # Convert degrees to radians
        inc = np.deg2rad(inc)
        az = np.deg2rad(90+az)  # convert b/c measured relative to north

        # Incidence-based components
        self.vFactor = np.cos(inc)  # vertical
        self.hFactor = np.sin(inc)  # horizontal

        # Horizontal components
        self.eFactor = np.cos(az)  # east
        self.nFactor = np.sin(az)  # north


def isce_geom_to_vector(inc, az, verbose=False):
    '''
    Convert the incidence and azimuth angles to vectors based on the ISCE
     convention and scaling factors given by the isce_sensitivity_factors class.

    INPUTS
        inc and az are from the ISCE los.rdr.geo file, and are given in degrees
    '''
    # Determine scaling factors
    S = isce_sensitivity_factors(inc, az)

    # Formulate pointing vector
    Px = S.hFactor*S.eFactor  # east component
    Py = S.hFactor*S.nFactor  # north component
    Pz = S.vFactor  # vertical component

    # Report if requested
    if verbose == True:
        print('Median pointing vector:')
        print('\tPx: {:.6f}'.format(np.median(Px)))
        print('\tPy: {:.6f}'.format(np.median(Py)))
        print('\tPz: {:.6f}'.format(np.median(Pz)))

    return Px, Py, Pz


def isce_to_los(ux, uy, uz, inc, az, verbose=False):
    '''
    Project 3D displacements into sensor line of sight based on the ISCE 
     convention and scaling factors given by the isce_sensitivity_factors class.
    '''
    # Determine scaling factors
    S = isce_sensitivity_factors(inc, az)

    # Compute line of sight displacement
    LOS = ux*S.eFactor + uy*S.nFactor + uz*S.vFactor

    # Report if requested
    if verbose == True:
        print('Median LOS value: {:.6f}'.format(np.median(LOS)))

    return LOS



### ISCE LOADING ---
def load_isce_los_dataset(fname, verbose=False):
    '''
    Load ISCE los.rdr.geo file as a GDAL data set, and retrieve incidence
     and azimuth maps.

    Inherits load_gdal_dataset from IOsupport and parse_isce_los.
    '''
    # Load data set
    LOS = load_gdal_dataset(fname, verbose=verbose)

    # Parse data set
    inc, az = parse_isce_los(LOS, verbose=verbose)

    return inc, az


def parse_isce_los(ISCEdataset, verbose=False):
    '''
    Correctly parse the incidence and azimuth data from an ISCE los.rdr.geo file.

    From the los.rdr.geo.xml file:
    Two channel Line-Of-Sight geometry image (all angles in degrees).
    Represents vector drawn from target to platform.
    Channel 1: Incidence angle measured from vertical at target (always +ve).
    Channel 2: Azimuth angle measured from North in Anti-clockwise direction.
    '''
    # Retrieve maps
    inc = ISCEdataset.GetRasterBand(1).ReadAsArray()
    az = ISCEdataset.GetRasterBand(2).ReadAsArray()

    return inc, az



### RECOMPOSITION ---
def invert_from_los(los, Px, Py, Pz, mask=None, nComponents=2, verbose=False):
    '''
    "Recompose" an observed signal (e.g., displacement or velocity) based on
     two or more lines of sight, following Wright et al. (2004).

    INPUTS
        los is a SxMxN array of LOS observations, where M and N are the
         north-south and east-west dimensions, and S is the number of
         obersvations (2-4).
        inc is a SxMxN array of incidence values
        az is a SxMxN array of azimuth values

    OUTPUTS
        E, (N,) Z are the inferred displacements/velocities, returned as
         MxN arrays

    First, convert the incidence and azimuth angles into look vectors using
    For each pixel, format a 2- or 3-component matrix, then invert for the 
     components of motion at that pixel.

    I.e., 
        P.u = r => u = (PT.P)-1 PT r

    where P is the matrix of pointing vectors, u is the true displacement
     vector, and r is the vector of observed range displacements

    Equivalently, for a weighted scheme, 
        W.P.u = W.r => u = (PT.W.P)-1 PT W r
    '''
    # Check all dimensions are correct
    nObs, M, N = inversion_checks(los, Px, Py, Pz, nComponents)

    if verbose == True:
        print('Solving for {:d} components based on {:d} observations'.format(nComponents, nObs))

    # Create matrices based on number of components for which to solve
    #  Pmatrices will be (MN x nObs x nComponents) in size
    if verbose == True: print('creating matrices')

    if nComponents == 2:
        Pmatrices = create_2comp_pointing_matrices(Px, Pz)

    elif nComponents == 3:
        Pmatrices = create_3comp_pointing_matrices(Px, Py, Pz)

    # Format los measurements into a (MN x nObs) array
    los = los.reshape(nObs, M*N)

    # Invert for displacement
    if verbose == True: print('solving for components')

    uhat = solve_los_for_components(Pmatrices, los)
    uhat = uhat.reshape(nComponents, M, N)  # reshape array

    return uhat


def inversion_checks(phs, Px, Py, Pz, nComponents):
    '''
    Check that inversion inputs are consistent.
    '''
    # Setup
    nObs, M, N = phs.shape  # number of observations

    # Run checks
    assert nObs <= 8, '{:d} observations detected'.format(nObs)
    assert phs.shape == Px.shape == Py.shape == Pz.shape, \
        '''All maps not equal in size:
Imgs {:d}x{:d}x{:d}
Px {:d}x{:d}x{:d}
Py {:d}x{:d}x{:d}
Pz {:d}x{:d}x{:d}'''.\
        format(*phs.shape, *Px.shape, *Py.shape, *Pz.shape)
    if nComponents ==3: assert nObs >= 3, \
        '{:d} oberservations is not sufficient to solve for three components of displacement'.\
        format(nObs)

    return nObs, M, N


def create_2comp_pointing_matrices(Px, Pz):
    '''
    Create 2-component pointing vector matrices, assuming zero sensitivity
     in the N-S direction.

    Matrices have the form:
        | p1x p1z |
    P = |  .   .  | 
        | pSx pSz |

    INPUTS
        Px, Pz are SxMxN arrays, where P is the number of observations and
         M and N are the map height and width.

    OUTPUTS
        matrices is a MxN-length list of matrices
    '''
    # Setup
    nObs, M, N = Px.shape
    nPixels = M*N
    Pmatrices = np.zeros((nObs, 2, nPixels))

    # Create matrices
    k = 0  # start counter
    for i in range(M):
        for j in range(N):
            Pmatrices[:,:,k] = np.column_stack([Px[:,i,j], Pz[:,i,j]])
            k += 1  # update counter

    return Pmatrices


def create_3comp_pointing_matrices(Px, Py, Pz):
    '''
    Create 3-component pointing vector matrices, assuming zero sensitivity
     in the N-S direction.

    Matrices have the form:
        | p1x p1y p1z |
    P = |  .   .   .  | 
        | pSx pSy pSz |

    INPUTS
        Px, Pz are SxMxN arrays, where P is the number of observations and
         M and N are the map height and width.

    OUTPUTS
        matrices is a MxN-length list of matrices
    '''
    # Setup
    nObs, M, N = Px.shape
    nPixels = M*N
    Pmatrices = np.zeros((nObs, 3, nPixels))

    # Create matrices
    k = 0  # start counter
    for i in range(M):
        for j in range(N):
            Pmatrices[:,:,k] = np.column_stack([Px[:,i,j], Py[:,i,j], Pz[:,i,j]])
            k += 1  # update counter

    return Pmatrices


def solve_los_for_components(Pmatrices, los):
    '''
    Estimate displacement given the pointing vectors in Pmatrices, and the LOS
     displacement vectors in los.
    
    INPUTS
        Pmatrices is a (nObs x nComponents x nPixels) array
        los is a (nComponents x nPixels) array
    '''
    # Setup
    assert Pmatrices.shape[2] == los.shape[1], \
        'Nb matrices ({:d}) must equal nb LOS ({:d})'.\
        format(Pmatrices.shape[2], los.shape[1])
    nPixels = Pmatrices.shape[2]  # how many data points
    nComponents = Pmatrices.shape[1]  # number of components

    # Empty displacements array
    uhat = np.zeros((nComponents, nPixels))

    # Solve each matrix
    for i in range(nPixels):
        try:
            # Retrieve pointing vector
            P = Pmatrices[:,:,i]

            # Retrieve LOS vector
            l = los[:,i]

            # Invert for solution
            uhat[:,i] = np.linalg.inv(np.dot(P.T, P)).dot(P.T).dot(l)
        except:
            print(P, 'Warning: Matrix is singular. Filling with zeros.')
            uhat[:,i] = np.zeros(nComponents)

    return uhat