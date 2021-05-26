#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Project GPS velocities or displacements into LOS.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import confirm_outdir, confirm_outname_ext, append_fname, load_gdal_dataset, GPSdata
from GeoFormatting import DS_to_extent, lola_to_xy
from RasterResampling import sample_points_from_raster
from RadarGeometries import aria_to_los, isce_to_los
from Masking import create_mask


### PARSER ---
Description = '''Project GPS velocities or displacements into satellite line of sight (LOS).
These routines work for ARIA and ISCE conventions.'''

Examples = ''''''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    GPSargs = parser.add_argument_group('GPS INPUTS')
    GPSargs.add_argument(dest='GPSname', type=str,
        help='GPS file name.')
    GPSargs.add_argument('-lon','--lon-column', dest='lonCol', type=int, default=0,
        help='Longitude column')
    GPSargs.add_argument('-lat','--lat-column', dest='latCol', type=int, default=1,
        help='Latitude column')
    GPSargs.add_argument('-e','--east-column', dest='Ecol', type=int, required=True,
        help='Column for east component (starts at 0, [None]).')
    GPSargs.add_argument('-n','--north-column', dest='Ncol', type=int, required=True,
        help='Column for north component (starts at 0, [None]).')
    GPSargs.add_argument('-z','--vert-column', dest='Vcol', type=int, default=None,
        help='Column for vertical component (starts at 0, [None]).')
    GPSargs.add_argument('--header-rows', dest='headerRows', type=int, default=0,
        help='Number of header rows to skip ([0], 1, 2, ...)')
    GPSargs.add_argument('--delimiter', dest='delimiter', type=str, default=' ',
        help='Text file delimiter ([\' \'],... etc.)')
    GPSargs.add_argument('-b','--bbox', dest='bbox', type=float, nargs=4, default=[-180, 180, -90, 90],
        help='Bounding box (WESN).')


    ProjectionArgs = parser.add_argument_group('PROJECTION INPUTS')
    ProjectionArgs.add_argument('-c','--convention', dest='convention', required=True,
        help='Convention (ARIA, ISCE).')
    ProjectionArgs.add_argument('-i','--incidence', dest='incInpt', type=str, default=None,
        help='Incidence (float or filename).')
    ProjectionArgs.add_argument('-a','--azimuth', dest='azInpt', type=str, default=None,
        help='Azimuth (float or filename).')
    ProjectionArgs.add_argument('-g','--geometry', dest='geomFile', type=str, default=None,
        help='ISCE geometry file.')
    ProjectionArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')
    ProjectionArgs.add_argument('--normalize-to-incidence', dest='norm2inc', action='store_true',
        help='Create a second file in which the GPS data are normalized by the incidence angle.')
    ProjectionArgs.add_argument('--scale-factor', dest='scaleFactor', type=float, default=1,
        help='Scale factor for GPS, e.g., 0.001 converts mm to m. ([1]).')


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



### GPS TO LOS ---
class GPStoLOS:
    def __init__(self, GPS, convention, incInpt=None, azInpt=None, geomFile=None, maskArgs=None, scaleFactor=1, verbose=False):
        '''
        Project GPS station measurements into satellite LOS.
        '''
        # Parameters
        self.convention = convention.lower()
        self.verbose = verbose

        # Check inputs are consistent
        self.__check_gps_inputs__(GPS)
        self.__check_geom_inputs__(incInpt, azInpt, geomFile)

        # Determine input types
        self.__determine_geom_types__(incInpt, azInpt, geomFile)

        if self.verbose == True: print('Initial checks passed')

        # Load values
        self.__load_geometries__(incInpt, azInpt, geomFile)

        # Check file sizes
        self.__check_file_sizes__()

        # Mask
        self.__compute_mask__(maskArgs)

        # Convert geographic coords to pixel values
        self.__compute_los__()

        # Scale by factor
        self.__scale_by_factor__(scaleFactor)

        # Compute LOS values
        if self.verbose == True: print('Projecting GPS into satellite line of sight (LOS)')


    ## Checks
    def __check_gps_inputs__(self, GPS):
        '''
        Check that GPS object has all the necessary attributes.
        If vertical component is not specified, 
        '''
        if self.verbose == True: print('Checking GPS inputs')

        # Check that vertical data are included, or assume zero values
        if GPS.V is None:
            if self.verbose == True: print('Vertical values not detected. Assuming zero displacement.')
            GPS.assign_zero_vertical()

        # Ascribe to object
        self.GPS = GPS


    def __check_geom_inputs__(self, incInpt, azInpt, geomFile):
        '''
        Check that file types are consistent.
        '''
        if self.verbose == True:
            print('Checking inputs are consistent with {:s} convention'.format(self.convention))

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

        # Report if requested
        if self.verbose == True: print('... input sizes consistent.')


    ## Inputs
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

        # Check if input is already a numpy array
        elif type(value) == np.ndarray:
            valueType = 'array'

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


    ## Masking
    def __compute_mask__(self, maskArgs):
        '''
        Compute a common mask for the geometries and inputs, if maps are specified.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Creating mask')

        # Create mask if necessary
        if self.geomType == 'map':
            # Initial mask
            self.mask = create_mask(self.inc, maskArgs, verbose=self.verbose)

        else:
            self.mask = None


    ## LOS projection
    def __format_geom_inputs__(self):
        '''
        Convert GPS coordinates to image pixel values.
        Crop and mask.
        '''
        if self.verbose == True: print('Converting GPS coordinates to pixel values')

        # Format inputs based on geometry type
        if self.geomType == 'float':
            # Use single values for incidence and azimuth
            self.incValues = self.inc
            self.azValues = self.az

        elif self.geomType in ['map', 'array']:
            # Sample from map values
            if self.verbose == True: print('Retreiving sample point values')

            # Crop to spatial extent of masks
            self.GPS.crop_to_extent(*self.extent)

            # Convert coordinates to pixel values
            px, py = lola_to_xy(self.tnsf, self.GPS.lon, self.GPS.lat, verbose=self.verbose)

            # Retreive sample points from raster
            sX, sY, self.incValues, sNdx = sample_points_from_raster(self.inc, px, py, mask=self.mask, verbose=self.verbose)
            sX, sY, self.azValues, sNdx = sample_points_from_raster(self.az, px, py, mask=self.mask, verbose=self.verbose)

            # Retain only non-masked GPS values
            self.GPS.crop_by_index(sNdx)


    def __compute_los__(self):
        '''
        Project the GPS station measurements into LOS.
        If geometry is given by single specified values, proceed directly.
        If geometry is provided as maps, sample those maps at the station locations.

        Wrapper for __compute_los_aria__ and __compute_los_isce__.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Projecting into LOS using {:s} convention'.format(self.convention))

        # Format geometry inputs
        self.__format_geom_inputs__()

        # Compute LOS
        if self.convention in ['aria']:
            self.los = aria_to_los(self.GPS.E, self.GPS.N, self.GPS.V, self.incValues, self.azValues, 
                verbose=self.verbose)

        elif self.convention in ['isce']:
            self.los = isce_to_los(self.GPS.E, self.GPS.N, self.GPS.V, self.incValues, self.azValues,
                verbose=self.verbose)

        # Number of stations projected
        self.nStations = len(self.los)


    def __scale_by_factor__(self, scaleFactor=1):
        '''
        Multiply by the given scale factor, e.g.,
            1 produces no change in scale
            0.001 converts mm to m
        '''
        if (self.verbose == True) and (scaleFactor != 1):
            print('*'*32)
            print('Scaling by factor: {:f}'.format(scaleFactor))

        # Scale LOS values
        self.los *= scaleFactor


    ## Normalization to sine of incidence
    def normalize_sine_incidence(self):
        '''
        Normalize the data to the sine of the incidence angle.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Normalizing by the sine of the incidence angle')

        # Normalize by sin(incidence)
        self.losNormd = self.los/np.sin(np.deg2rad(self.incValues))


    ## Saving
    def save(self, outName):
        '''
        Save the LOS measurements to a text file.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Saving to file')

        # Format output name
        outName = confirm_outname_ext(outName, ext=['txt'], verbose=self.verbose)
        confirm_outdir(outName)

        # Setup
        header = 'lon lat los\n'
        stationFmt = '{:f} {:f} {:f}\n'

        # Write to file
        with open(outName, 'w') as outFile:
            # Write header
            outFile.write(header)

            # Write stations to file
            [outFile.write(stationFmt.format(self.GPS.lon[i], self.GPS.lat[i], self.los[i])) \
                for i in range(self.nStations)]

        # Save sin(inc)-normalized data if available
        if hasattr(self, 'losNormd'):
            if self.verbose == True: print('Saving incidence-normalized data')

            outName = append_fname(outName, '_inc-normd', verbose=self.verbose)
            with open(outName, 'w') as outFile:
                # Write header
                outFile.write(header)

                # Write stations to file
                [outFile.write(stationFmt.format(self.GPS.lon[i], self.GPS.lat[i], self.losNormd[i])) \
                for i in range(self.nStations)]


    ## Plotting
    def plot(self):
        '''
        Plot the GPS stations inputs and outputs.
        '''
        if self.verbose == True: print('Plotting inputs and outputs')

        # Spawn figure and axes
        fig, [[axE, axN], [axV, axLOS]] = plt.subplots(figsize=(8,8), nrows=2, ncols=2)

        # Determine color scale bounds
        vmin, vmax = self.__gps_colorscale__()

        # Plot east measurements
        self.__plot_gps__(fig, axE, self.GPS.E, title='East', vmin=vmin, vmax=vmax)
        self.__plot_gps__(fig, axN, self.GPS.N, title='North', vmin=vmin, vmax=vmax)
        self.__plot_gps__(fig, axV, self.GPS.V, title='Vertical', vmin=vmin, vmax=vmax)
        self.__plot_gps__(fig, axLOS, self.los, title='LOS', vmin=None, vmax=None)


    def __gps_colorscale__(self):
        '''
        Determine the min/max color bounds to use on the GPS data.
        Clip to +- 1 percent.
        '''
        # Concatenate GPS values
        data = np.concatenate((self.GPS.E, self.GPS.N, self.GPS.V, self.los), axis=0)

        # Find 1st and 99th percentiles
        vmin, vmax = np.percentile(data, [1, 99])

        return vmin, vmax


    def __plot_gps__(self, fig, ax, data, title=None, vmin=None, vmax=None):
        '''
        Subplot for GPS measurements.
        '''
        # Plot data
        cax = ax.scatter(self.GPS.lon, self.GPS.lat, 9, c=data,
            cmap='viridis', vmin=vmin, vmax=vmax)

        # Format subplot
        ax.set_title(title)
        ax.set_aspect(1)
        fig.colorbar(cax, ax=ax, orientation='horizontal')



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Load data
    # Instantiate GPS object
    GPS = GPSdata(verbose=inps.verbose)

    # Load GPS data
    GPS.load_from_file(fname=inps.GPSname,
        lonCol=inps.lonCol, latCol=inps.latCol,
        Ecol=inps.Ecol, Ncol=inps.Ncol, Vcol=inps.Vcol,
        headerRows=inps.headerRows, delimiter=inps.delimiter)

    # Crop to extent
    GPS.crop_to_extent(*inps.bbox)


    ## LOS projection
    LOS = GPStoLOS(GPS, inps.convention,
        incInpt=inps.incInpt, azInpt=inps.azInpt, geomFile=inps.geomFile,
        maskArgs=inps.maskArgs,
        scaleFactor=inps.scaleFactor,
        verbose=inps.verbose)

    # Normalize to incidence if requested
    if inps.norm2inc == True:
        LOS.normalize_sine_incidence()


    ## Outputs
    # Save to file
    LOS.save(inps.outName)

    # Plot if requested
    if inps.plot == True:
        LOS.plot()


    plt.show()