#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Extract the three-component ground-to-sensor look vectors for ARIA or ISCE data.

FUTURE IMPROVEMENTS

TESTING STATUS
In development
'''

### IMPORT MODULES ---
import argparse
import os
import matplotlib.pyplot as plt
from IOsupport import load_gdal_dataset, confirm_outdir, confirm_outname_ext, save_gdal_dataset
from Checks import check_dataset_sizes
from GeoFormatting import DS_to_extent
from Masking import create_mask
from RadarGeometries import aria_geom_to_vector, isce_geom_to_vector
from Viewing import raster_multiplot


### PARSER ---
Description = '''Extract the three-component ground-to-sensor look vectors for ARIA or ISCE data.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')

    InputArgs.add_argument('-c','--convention', dest='convention', required=True,
        help='Convention (ARIA, ISCE).')

    InputArgs.add_argument('-i','--incidence', dest='incFile', type=str, default=None,
        help='ARIA incidence files.')
    InputArgs.add_argument('-a','--azimuth', dest='azFile', type=str, default=None,
        help='ARIA azimuth files.')
    InputArgs.add_argument('-g','--geometry', dest='geomFile', type=str, default=None,
        help='ISCE geometry files.')

    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true', 
        help='Plot inputs and outputs.')
    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### LOOK VECTORS CLASS ---
class look_extraction:
    def __init__(self, convention, incFile, azFile, geomFile, verbose=False):
        '''
        Extract look vectors from geometry, azimuth, and/or incidence files of
         ISCE or ARIA data.

        Inherits
            os
            IOsupport: confirm_outdir, confirm_outname_ext,
             load_gdal_dataset save_gdal_dataset
            Checks: check_dataset_sizes
            GeoFormatting: DS_to_extent
            Masking: create_mask
            RadarGeometries: aria_geom_to_vector, isce_geom_to_vector
            Viewing: raster_multiplot

        To initialize, check that inputs are consistent.
        '''
        # Parameters
        self.convention = convention.lower()
        self.verbose = verbose

        # Check that files are valid and suitable for the convention
        self.__check_inputs__(azFile, incFile, geomFile)

        # Load geometry data
        self.__load_inputs__(azFile, incFile, geomFile)

        # Create mask
        self.__create_mask__()

        # Convert geometry angles to vectors
        self.__convert_to_vectors__()


    ## Checks
    def __check_inputs__(self, incFile, azFile, geomFile):
        '''
        Check that files are valid and suitable for the specified convention.
        Wrapper for __check_aria_inputs__ and __check_isce_inputs__.
        '''
        if self.verbose == True: print('Checking inputs for {:s} convention'.format(self.convention.upper()))

        if self.convention in ['aria']:
            # Check ARIA inputs
            self.__check_aria_inputs__(incFile, azFile)

        elif self.convention in ['isce']:
            # Check ISCE inputs
            self.__check_isce_inputs__(geomFile)


    def __check_aria_inputs__(self, incFile, azFile):
        '''
        Check inputs based on ARIA conventions.
        '''
        # Check that appropriate files are specified
        assert azFile is not None and incFile is not None, \
            'Azimuth file and incidence file must be specified for ARIA data.'

        # Check that specified files are valid
        assert os.path.exists(incFile), 'Incidence file does not exist!'
        assert os.path.exists(azFile), 'Azimuth file does not exist!'


    def __check_isce_inputs__(self, geomFile):
        '''
        Check inputs based on ISCE conventions.
        '''
        # Check that appropriate files are specified
        assert geomFile is not None, \
            'Geometry file must be specified for ISCE data.'

        # Check that specified file is valid
        assert os.path.exists(geomFile), 'Geometry file does not exist!'


    ## Loading
    def __load_inputs__(self, incFile, azFile, geomFile):
        '''
        Load the input geometries.
        Wrapper for __load_aria_geom__ and __load_isce_geom__.
        '''
        if self.verbose == True: print('Loading {:s} data'.format(self.convention.upper()))

        if self.convention in ['aria']:
            # Load ARIA data
            self.__load_aria_geom__(incFile, azFile)

        elif self.convention in ['isce']:
            # Load ISCE data
            self.__load_isce_geom__(geomFile)


    def __load_aria_geom__(self, incFile, azFile):
        '''
        Load the incidence and azimuth data from ARIA files.
        '''
        # Load incidence
        DSinc = load_gdal_dataset(incFile, verbose=self.verbose)
        self.inc = DSinc.GetRasterBand(1).ReadAsArray()

        # Load azimuth
        DSaz = load_gdal_dataset(azFile, verbose=self.verbose)
        self.az = DSaz.GetRasterBand(1).ReadAsArray()

        # Check dataset sizes
        check_dataset_sizes([DSinc, DSaz])

        # Parse spatial data
        self.tnsf = DSinc.GetGeoTransform()
        self.proj = DSinc.GetProjection()
        self.extent = DS_to_extent(DSinc, verbose=self.verbose)


    def __load_isce_geom__(self, geomFile):
        '''
        Load the incidence and azimuth data from ISCE LOS file.
        '''
        # Load multi-band geometry file
        DSgeom = load_gdal_dataset(geomFile, verbose=self.verbose)
        self.inc = DSgeom.GetRasterBand(1).ReadAsArray()
        self.az = DSgeom.GetRasterBand(2).ReadAsArray()

        # Parse spatial data
        self.tnsf = DSgeom.GetGeoTransform()
        self.proj = DSgeom.GetProjection()
        self.extent = DS_to_extent(DSgeom, verbose=self.verbose)


    ## Masking
    def __create_mask__(self):
        '''
        Create masking using the Masking: create_mask function.
        '''
        self.mask = create_mask(self.inc, ['bg'])


    ## Vector conversion
    def __convert_to_vectors__(self):
        '''
        Convert incidence and azimuth angles to vectors based on the specified
         convention.
        Wrapper for __geom_to_vectors_aria__ and __geom_to_vectors_isce__.
        '''
        if self.verbose == True: print('Converting {:s} geometry to vectors'.format(self.convention.upper()))

        if self.convention in ['aria']:
            # Check ARIA inputs
            self.__geom_to_vectors_aria__()

        elif self.convention in ['isce']:
            # Check ISCE inputs
            self.__geom_to_vectors_isce__()


    def __geom_to_vectors_aria__(self):
        '''
        Convert angles to vectors following the ARIA convention.
        '''
        self.px, self.py, self.pz = aria_geom_to_vector(self.inc, self.az, verbose=self.verbose)


    def __geom_to_vectors_isce__(self):
        '''
        Convert angles to vectors following the ISCE convention.
        '''
        self.px, self.py, self.pz = isce_geom_to_vector(self.inc, self.az, verbose=self.verbose)


    ## Saving
    def save(self, outName):
        '''
        Save the vector components fields to a three-band GeoTiff using GDAL.
        '''
        if self.verbose == True: print('Saving vector fields')

        # Check outname
        confirm_outdir(outName)  # confirm output directory exists
        outName = confirm_outname_ext(outName)  # confirm output extension if GeoTiff

        # Save to GDAL data set
        save_gdal_dataset(outName, [self.px, self.py, self.pz], mask=self.mask,
            proj=self.proj, tnsf=self.tnsf, verbose=self.verbose)


    ## Plotting
    def plot(self):
        '''
        Plot input geometry maps and output LOS vector fields.
        '''
        # Plot incidence and azimuth data
        raster_multiplot([self.inc, self.az], ncols=2,
            mask=self.mask, extent=self.extent,
            cbarOrient='horizontal',
            titles=['Incidence', 'Azimuth'])

        # Plot vector components
        raster_multiplot([self.px, self.py, self.pz], ncols=3,
            mask=self.mask, extent=self.extent,
            cbarOrient='horizontal',
            titles=['Pointing x', 'Pointing y', 'Pointing z'])



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Compute look vectors
    # Compute vector fields
    looks = look_extraction(convention=inps.convention,
        incFile=inps.incFile, azFile=inps.azFile, geomFile=inps.geomFile,
        verbose=inps.verbose)

    # Save to file
    looks.save(inps.outName)

    # Plot if specified
    if inps.plot == True: looks.plot()


    plt.show()