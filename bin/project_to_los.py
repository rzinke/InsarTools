#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Project three components of a field to line of sight.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import load_gdal_dataset, confirm_outdir, confirm_outname_ext, save_gdal_dataset
from GeoFormatting import DS_to_extent
from RadarGeometries import aria_to_los, aria_geom_to_vector, isce_to_los, isce_geom_to_vector
from Masking import create_mask
from Viewing import plot_look_vectors, raster_multiplot


### PARSER ---
Description = '''Project three components of motion into satellite line of sight (LOS).
These routines work for ARIA and ISCE conventions.'''

Examples = '''# Provide ARIA look direction as two values, and motion as three values
project_to_los.py
'''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')

    InputArgs.add_argument('-c','--convention', dest='convention', required=True,
        help='Convention (ARIA, ISCE).')

    InputArgs.add_argument('-a','--azimuth', dest='azInpt', type=str, default=None,
        help='Azimuth (float or filename).')
    InputArgs.add_argument('-i','--incidence', dest='incInpt', type=str, default=None,
        help='Incidence (float or filename).')
    InputArgs.add_argument('-g','--geometry', dest='geomFile', type=str, default=None,
        help='ISCE geometry file.')
    InputArgs.add_argument('-e','--east', dest='Einpt', type=str, required=True,
        help='East component (float or filename).')
    InputArgs.add_argument('-n','--north', dest='Ninpt', type=str, required=True,
        help='North component (float or filename).')
    InputArgs.add_argument('-z','--vert', dest='Vinpt', type=str, required=True,
        help='Vertical component (float or filename).')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')

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



### CHECKS ---
class LOSproject:
    def __init__(self, convention, incInpt, azInpt, geomFile, Einpt, Ninpt, Vinpt, maskArgs, verbose=False):
        '''
        Class for converting three components of motion to satellite line of 
         sight based on radar geometry.

        Inherits:
            os
            numpy
            IOsupport: load_gdal_dataset, confirm_outdir, confirm_outname_ext, save_gdal_dataset
            GeoFormatting: DS_to_extent
            RadarGeometries: aria_to_los, aria_geom_to_vector, isce_to_los, isce_geom_to_vector
            Masking create_mask
            Viewing: plot_look_vectors, raster_multiplot

        To initialize, check that inputs types are consistent. Then compute LOS.
        '''
        # Parameters
        self.convention = convention.lower()
        self.verbose = verbose

        # Check inputs are consistent with convention
        self.__check_inputs__(incInpt, azInpt, geomFile)

        # Determine input types
        self.__determine_geom_types__(incInpt, azInpt, geomFile)
        self.__determine_displacement_types__(Einpt, Ninpt, Vinpt)

        if self.verbose == True: print('Initial checks passed')

        # Load values
        self.__load_geometries__(incInpt, azInpt, geomFile)
        self.__load_displacements__(Einpt, Ninpt, Vinpt)

        # Check file sizes
        self.__check_file_sizes__()

        # Mask
        self.__compute_mask__(maskArgs)

        # Compute LOS values
        self.__compute_los__()


    ## Inputs
    def __check_inputs__(self, incInpt, azInpt, geomFile):
        '''
        Check that given inputs are consistent with convention.
        '''
        if self.verbose == True:
            print('Checking inputs are consistent with {:s} convention'.format(self.convention.upper()))

        # Check convention is valid
        assert self.convention in ['aria', 'isce'], 'Convention must be ARIA or ISCE'

        # For ARIA convention
        if self.convention in ['aria']:
            assert incInpt is not None and azInpt is not None, \
                'Incidence and azimuth angles/files must be specified for ARIA convention'

        # For ISCE convention
        elif self.convention in ['isce']:
            # Void incidence and azimuth if geometry file is specified
            if geomFile is not None:
                if incInpt is not None or azInpt is not None:
                    print('ISCE geometry file specified, overriding incidence and azimuth values')
                incInpt = None
                azInpt = None
            else:
                assert incInpt is not None and azInpt is not None, \
                    'ISCE los.rdr.geo file must be given, otherwise incidence and azimuth angles must be specified in ISCE geometries'


    def __determine_geom_types__(self, incInpt, azInpt, geomFile):
        '''
        Determine the types of files specified and check they are consistent
         with the convention and with each other.

        Note geomType is a generic term for all geometry types.
        '''
        if self.verbose == True: print('Checking geometry inputs are consistent type')

        # Check incidence value
        if incInpt is not None:
            incType = self.__determine_type__(incInpt, label='incidence')

        if azInpt is not None:
            azType = self.__determine_type__(azInpt, label='azimuth')

            # Ensure types are the same
            assert incType == azType, 'Incidence and azimuth input types must be the same (float/map)'

            # Generic type for geometry
            self.geomType = incType

        if geomFile is not None:
            self.geomType = self.__determine_type__(geomFile, label='geometry')

            # Ensure proper formatting for geometry file
            assert self.geomType == 'map', 'Geometry file must be GDAL-readable data set'


    def __determine_displacement_types__(self, Einpt, Ninpt, Vinpt):
        '''
        Determine the types of files specified and check they are consistent
         with each other.
        '''
        if self.verbose == True: print('Checking displacement inputs are consistent type')

        # Check east value
        Etype = self.__determine_type__(Einpt, label='east')

        # Check north value
        Ntype = self.__determine_type__(Ninpt, label='north')

        # Check vertical value
        Vtype = self.__determine_type__(Vinpt, label='vertical')

        # Check that input types are the same
        assert Etype == Ntype == Vtype, 'Displacement input types must be the same (float/map)'

        # Generic type for displacement
        self.inptType = Etype


    def __determine_type__(self, value, label=None):
        '''
        Determine the input type (float or map) for the given input.
        '''
        # Test if value is an available file
        if os.path.exists(value):
            # Treat as file and load file using GDAL
            valueType = 'map'

            if self.verbose == True and label is not None:
                print('{:s} value identified as a map file'.format(label))

        # Test if value type can be converted to float
        else:
            try:
                # Try using float value
                float(value)
                valueType = 'float'

                if self.verbose == True and label is not None:
                    print('{:s} value identified as float'.format(label))

            except:
                if label is not None:
                    print('{:s} value could not be identified. Exiting.'.format(label))
                else:
                    print('Value could not be identified. Exiting.')

                exit()

        return valueType


    ## Loading
    def __load_geometries__(self, incInpt, azInpt, geomFile):
        '''
        Load incidence and azimuth angles based on the convention and value
         types.
        Wrapper for __load_aria_geoms__ and __load_isce_geoms__.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Loading {:s} geometries'.format(self.convention.upper()))

        # Load ARIA geometries
        if self.convention in ['aria']:
            self.__load_aria_geoms__(incInpt, azInpt)

        # Load ISCE geometries
        elif self.convention in ['isce']:
            self.__load_isce_geoms__(incInpt, azInpt, geomFile)


    def __load_aria_geoms__(self, incInpt, azInpt):
        '''
        Load or format ARIA input data.
        ARIA incidence and azimuth data are always specified separately, so this
         funciton effectively acts as a wrapper for __load_inc_az__.
        '''
        # Retrieve incidence and azimuth values
        self.inc, self.az = self.__load_inc_az__(incInpt, azInpt)


    def __load_isce_geoms__(self, incInpt, azInpt, geomFile):
        '''
        Load or format ISCE data.
        If separate incidence and azimuth values are specified, retrieve using
         the __load_inc_az__ subroutine. Otherwise, load and parse the ISCE
         geometry file.
        '''
        # Format or load incidence and azimuth
        if incInpt is not None and azInpt is not None:
            # Load values specified separately
            self.inc, self.az = self.__load_inc_az__(incInpt, azInpt)

        elif geomFile is not None:
            # Load ISCE geometry data set
            DSgeom = load_gdal_dataset(geomFile, verbose=self.verbose)
            self.inc = DSgeom.GetRasterBand(1).ReadAsArray()
            self.az = DSgeom.GetRasterBand(2).ReadAsArray()

            # Spatial parameters
            self.extent = DS_to_extent(DSgeom)
            self.proj = DSgeom.GetProjection()
            self.tnsf = DSgeom.GetGeoTransform()


    def __load_inc_az__(self, incInpt, azInpt):
        '''
        Load incidence and azimuth data based on file type.
        '''
        # Format or load incidence and azimuth
        if self.geomType == 'float':
            inc = float(incInpt)  # incidence value
            az = float(azInpt)  # azimuth value
        elif self.geomType == 'map':
            DSinc = load_gdal_dataset(incInpt, verbose=self.verbose)
            incExtent = DS_to_extent(DSinc)  # spatial extent
            inc = DSinc.GetRasterBand(1).ReadAsArray()  # map of incidence values

            DSaz = load_gdal_dataset(azInpt, verbose=self.verbose)
            azExtent = DS_to_extent(DSaz)  # spatial extent
            az = DSaz.GetRasterBand(1).ReadAsArray()  # map of azimuth values

            # Check spatial extents are the same
            assert incExtent == azExtent, 'Spatial extents must be the same'
            self.extent = incExtent  # commit to attributes
            self.proj = DSinc.GetProjection()
            self.tnsf = DSinc.GetGeoTransform()

        return inc, az


    def __load_displacements__(self, Einpt, Ninpt, Vinpt):
        '''
        Load displacement values based on the value types.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Loading displacements')

        # Load displacement values
        self.E = self.__load_displacement__(Einpt, self.inptType)
        self.N = self.__load_displacement__(Ninpt, self.inptType)
        self.V = self.__load_displacement__(Vinpt, self.inptType)


    def __load_displacement__(self, inpt, inptType):
        '''
        Format or load value for displacement field based on file type.
        '''
        # Load the input based on file type
        if inptType == 'float':
            value = float(inpt)
        elif inptType == 'map':
            DS = load_gdal_dataset(inpt, verbose=self.verbose)
            extent = DS_to_extent(DS)  # spatial extent
            value = DS.GetRasterBand(1).ReadAsArray()

            if hasattr(self, 'extent'):
                assert extent == self.extent, 'Spatial extents must be the same'
            else:
                self.extent = extent
                self.proj = DS.GetProjection()
                self.tnsf = DS.GetGeoTransform()

        return value


    ## Checks
    def __check_file_sizes__(self):
        '''
        Check that file sizes are consistent between all geometries and displacements.
        '''
        if self.verbose == True: print('Checking file sizes are consistent')

        # Check geometry sizes
        if self.geomType == 'map':
            incShape = self.inc.shape
            azShape = self.az.shape
            assert incShape == azShape, \
                'Incidence and azimuth files must be the same size ({:d} x {:d}) and ({:d} x {:d})'.\
                format(*incShape, *azShape)

        # Check input sizes
        if self.inptType == 'map':
            Eshape = self.E.shape
            Nshape = self.N.shape
            Vshape = self.V.shape
            assert Eshape == Nshape == Vshape, \
                'E, N, and V files must be the same size ({:d} x {:d}) and ({:d} x {:d}) and ({:d} x {:d})'.\
                format(*Eshape, *Nshape, *Vshape)

        # Check geometry sizes match input sizes
        if (self.geomType == 'map') and (self.inptType == 'map'):
            assert self.inc.shape == self.E.shape, \
            'Geometry and displacement input sizes must be the same size ({:d} x {:d}) and ({:d} x {:d})'.\
            format(*self.inc.shape, *self.E.shape)

        # Report if requested
        if self.verbose == True: print('... input sizes consistent.')


    ## Masking
    def __compute_mask__(self, maskArgs):
        '''
        Compute a common mask for the geometries and inputs, if maps are specified.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Creating mask')

        # Create mask if necessary
        if (self.geomType == 'map') or (self.inptType == 'map'):
            # Initial mask
            self.mask = np.ones(self.E.shape)  # pass all

            if self.geomType == 'map':
                # Create mask based on the geometry inputs
                self.mask[create_mask(self.inc, maskArgs) == 0] = 0

            elif self.inptType == 'map':
                # Create mask based on the displacement inputs
                self.mask[create_mask(self.E, maskArgs) == 0] = 0
        else:
            self.mask = None


    ## LOS computation
    def __compute_los__(self):
        '''
        Compute the LOS value(s) based on the satellite geometry.
        This function is a wrapper for __aria_to_los__ and __isce_to_los__.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Computing LOS using the {:s} convention'.format(self.convention.upper()))

        # Convert to LOS using the ARIA convention
        if self.convention in ['aria']:
            self.LOS = aria_to_los(self.E, self.N, self.V, self.inc, self.az, verbose=self.verbose)

        # Convert to LOS using the ISCE convention
        elif self.convention in ['isce']:
            self.LOS = isce_to_los(self.E, self.N, self.V, self.inc, self.az, verbose=self.verbose)


    ## Saving
    def save(self, outName):
        '''
        Save the vector components fields to a three-band GeoTiff using GDAL.
        '''
        # Save to file if a map was computed
        if type(self.LOS) == np.ndarray:
            if self.verbose == True:
                print('*'*32)
                print('Saving to file')

            # Check outname
            confirm_outdir(outName)  # confirm output directory exists
            outName = confirm_outname_ext(outName)  # confirm output extension if GeoTiff

            # Save to GDAL data set
            save_gdal_dataset(outName, self.LOS, mask=self.mask,
                proj=self.proj, tnsf=self.tnsf, verbose=self.verbose)

        else:
            print('LOS value: {:.8f}'.format(self.LOS))


    ## Plotting
    def plot(self):
        '''
        Plot the geometry, inputs, and outputs based on type.
        Wrapper for __plot_geometries__ and __plot_inputs__.
        '''
        # Plot geometry
        self.__plot_geometries__()

        # Plot inputs
        self.__plot_inputs__()


    def __plot_geometries__(self):
        '''
        Plot setup geometries based on input type.
        '''
        # Plot geometry
        if self.geomType == 'float':
            # Determine look vector
            if self.convention in ['aria']:
                Px, Py, Pz = aria_geom_to_vector(self.inc, self.az)
            elif self.convention in ['isce']:
                Px, Py, Pz = isce_geom_to_vector(self.inc, self.az)

            fig, axes = plot_look_vectors(Px, Py, Pz)
            fig.suptitle('Geometry')

        elif self.geomType == 'map':
            # Plot rasters
            raster_multiplot([self.inc, self.az], ncols=2,
                mask=self.mask, extent=self.extent,
                cbarOrient='horizontal',
                titles=['incidence', 'azimuth'], suptitle='Geometry')


    def __plot_inputs__(self):
        '''
        Plot displacements based on input type.
        '''
        # Plot components and LOS
        if self.inptType == 'map':
            # Plot rasters
            raster_multiplot([self.E, self.N, self.V, self.LOS], mrows=2, ncols=2,
                mask=self.mask, extent=self.extent,
                cbarOrient='vertical',
                titles=['east', 'north', 'vertical', 'LOS'])



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Project into LOS
    projection = LOSproject(convention=inps.convention,
        azInpt=inps.azInpt, incInpt=inps.incInpt, geomFile=inps.geomFile,
        Einpt=inps.Einpt, Ninpt=inps.Ninpt, Vinpt=inps.Vinpt,
        maskArgs=inps.maskArgs,
        verbose=inps.verbose)

    # Save
    projection.save(inps.outName)

    # Plot if requested
    if inps.plot == True:
        projection.plot()

        plt.show()