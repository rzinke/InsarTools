'''
SHORT DESCRIPTION
Conversion factors and routines for different radar geometry conventions.

FUTUTRE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import numpy as np
from IOsupport import load_gdal_dataset


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


def aria_geom_to_los(ux, uy, uz, inc, az, verbose=False):
    '''
    Project 3D displacements into sensor line of sight based on the ARIA 
     convention and scaling factors given by the aria_sensitivity_factors class.
    '''
    # Determine scaling factors
    S = aria_sensitivity_factors(inc, az)

    # Compute line of sight displacement
    LOS = ux*S.eFactor + uy*S.nFactor + uz*S.vFactor

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


def isce_geom_to_los(ux, uy, uz, inc, az, verbose=False):
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