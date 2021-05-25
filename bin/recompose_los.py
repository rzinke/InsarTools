#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Convert two or more lines of sight into movement components.

FUTURE IMPROVEMENTS
    * Weighting for inversion scheme

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import load_gdal_datasets, confirm_outdir, confirm_outname_ext, save_gdal_dataset
from RasterResampling import gdal_resample, match_rasters
from Checks import check_dataset_sizes
from Viewing import plot_raster, raster_multiplot
from GeoFormatting import DS_to_bounds, DS_to_extent
from Masking import create_mask
from RadarGeometries import aria_geom_to_vector, isce_geom_to_vector, invert_from_los, parse_isce_los


### PARSER ---
Description = '''Project three components of motion into satellite line of sight (LOS).
These routines work for ARIA and ISCE conventions.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')

    InputArgs.add_argument('-c','--convention', dest='convention', required=True,
        help='Convention (ARIA, ISCE).')

    InputArgs.add_argument('-f','--ifg', dest='imgFiles', type=str, nargs='+', required=True,
        help='LOS file names.')
    InputArgs.add_argument('-i','--incidence', dest='incFiles', type=str, nargs='+', default=None,
        help='ARIA incidence files.')
    InputArgs.add_argument('-a','--azimuth', dest='azFiles', type=str, nargs='+', default=None,
        help='ARIA azimuth files.')
    InputArgs.add_argument('-g','--geometry', dest='geomFiles', type=str, nargs='+', default=None,
        help='ISCE geometry files.')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')
    InputArgs.add_argument('-n','--n-components', dest='nComponents', type=int, default=2,
        help='Number of displacement components for which to solve ([2], 3).')
    InputArgs.add_argument('--max-cpus', dest='maxCPUs', type=int, default=8,
        help='Maximum number of CPUs to use (1, 2, ..., [8], ...).')

    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('--plot-inputs', dest='plotInputs', action='store_true', 
        help='Plot inputs.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true', 
        help='Plot outputs.')
    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### RECOMPOSITION CLASS ---
class LOSrecomposition:
    def __init__(self, convention, imgFiles, incFiles, azFiles, geomFiles, maskArgs=[], nComponents=2, maxCPUs=8,
        verbose=False):
        '''
        Use the satellite geometry files and LOS values to "recompose" the
         EW, (NS,) and vertical components of motion.
        '''
        # Parameters
        self.convention = convention.lower()
        self.nComponents = nComponents
        self.maxCPUs = maxCPUs
        self.verbose = verbose

        # Check that inputs are consistent
        self.__check_inputs__(imgFiles, azFiles, incFiles, geomFiles)

        # Load data from files
        self.__load_data__(imgFiles, azFiles, incFiles, geomFiles)

        # Mask
        self.__create_mask__(maskArgs)

        # Convert geometries to look vectors
        self.__geom_to_vectors__()

        # Calculate displacement
        self.__recompose__()


    ## Checks
    def __check_inputs__(self, imgFiles, incFiles, azFiles, geomFiles):
        '''
        Check that the number of files is consistent.
        Store attribute nObs.
        '''
        if self.verbose == True: print('Checking input consistency...')

        # Check number of input files
        if self.convention in ['isce']:
            # ... for ISCE data
            assert len(imgFiles) == len(geomFiles), \
                'Number of image files ({:d}) must equal number of geometry files ({:d}).'.\
                format(len(imgFiles), len(geomFiles))

        elif self.convention in ['aria']:
            # ... for ARIA data
            assert len(imgFiles) == len(incFiles) == len(azFiles), \
                'Number of image files ({:d}) must equal number of incidence ({:d}) and azimuth ({:d}) files.'.\
                format(len(imgFiles), len(incFiles), len(azFiles))

        # Number of files for quick reference
        self.nObs = len(imgFiles)

        # Report if requested
        if self.verbose == True: print('\tnumber of input files is consistent: ({:d})'.format(self.nObs))


    ## Load data
    def __load_data__(self, imgFiles, incFiles, azFiles, geomFiles):
        '''
        Load data from given files.
        This is a wrapper for __load_data_isce__ and __load_data_aria__.
        '''
        if self.verbose == True:
                print('*'*32)
                print('Loading {:s} data'.format(self.convention.upper()))

        # Load ISCE data
        if self.convention in ['isce']:
            self.__load_data_isce__(imgFiles, geomFiles)

        # Load ARIA data
        elif self.convention in ['aria']:
            self.__load_data_aria__(imgFiles, azFiles, incFiles)


    def __load_data_aria__(self, imgFiles, incFiles, azFiles):
        '''
        Load ARIA data as GDAL data sets.
        '''
        # Load LOS displacement/velocity data
        imgDatasets = load_gdal_datasets(imgFiles, verbose=self.verbose)

        # Load geometry data sets
        incDatasets = load_gdal_datasets(incFiles, verbose=self.verbose)
        azDatasets = load_gdal_datasets(azFiles, verbose=self.verbose)

        # Save data set names for posterity
        self.dsNames = list(imgDatasets.keys())

        # Resample phase data sets to same size
        imgDatasets = match_rasters(imgDatasets, cropping='intersection', verbose=self.verbose)
        Mimg, Nimg = check_dataset_sizes(imgDatasets)

        # Store spatial parameters
        bounds = DS_to_bounds(imgDatasets[self.dsNames[0]])
        self.extent = DS_to_extent(imgDatasets[self.dsNames[0]])
        self.tnsf = imgDatasets[self.dsNames[0]].GetGeoTransform()
        self.proj = imgDatasets[self.dsNames[0]].GetProjection()

        # Resample geometry data sets to match phase data sets
        incDatasets = [gdal_resample(incDatasets[dsName], bounds=bounds, M=Mimg, N=Nimg) for dsName in dsNames]
        azDatasets = [gdal_resample(azDatasets[dsName], bounds=bounds, M=Mimg, N=Nimg) for dsName in dsNames]

        # Ensure data sets are the same size
        Minc, Ninc = check_dataset_sizes(incDatasets)
        Maz, Naz = check_dataset_sizes(azDatasets)

        assert (Mimg, Nimg) == (Minc, Ninc) == (Maz, Naz), \
            'Image data set dimensions ({:d} x {:d}) must equal incidence ({:d} x {:d}) and azimuth ({:d} x {:d})'.\
            format(Mimg, Nimg, Minc, Ninc, Maz, Naz)

        # Retrieve phase information from IFG data
        phsMaps = [imgDatasets[dsName].GetRasterBand(1).ReadAsArray() for dsName in imgDatasets.keys()]

        # Retrieve incidence and azimuth angles
        incMaps = [incDatasets[dsName].GetRasterBand(1).ReadAsArray() for dsName in incDatasets.keys()]
        azMaps = [azDatasets[dsName].GetRasterBand(1).ReadAsArray() for dsName in azDatasets.keys()]

        # Convert to arrays
        self.phs = np.array(phsMaps)
        self.inc = np.array(incMaps)
        self.az = np.array(azMaps)


    def __load_data_isce__(self, imgFiles, geomFiles):
        '''
        Load ISCE data as GDAL data sets.
        '''
        # Load LOS displacement/velocity data
        imgDatasets = load_gdal_datasets(imgFiles, verbose=self.verbose)

        # Load geometry data sets
        geomDatasets = load_gdal_datasets(geomFiles, verbose=self.verbose)

        # Save data set names for posterity
        self.dsNames = list(imgDatasets.keys())

        # Resample data sets to same size
        imgDatasets = match_rasters(imgDatasets, cropping='intersection', verbose=self.verbose)
        geomDatasets = match_rasters(geomDatasets, cropping='intersection', verbose=self.verbose)

        # Ensure data sets are the same size
        Mimg, Nimg = check_dataset_sizes(imgDatasets)
        Mgeom, Ngeom = check_dataset_sizes(geomDatasets)

        assert (Mimg, Nimg) == (Mgeom, Ngeom), \
            'Image data set dimensions ({:d} x {:d}) must equal geomtry data set ({:d} x {:d})'.\
            format(Mimg, Nimg, Mgeom, Ngeom)

        # Store spatial parameters
        self.extent = DS_to_extent(imgDatasets[self.dsNames[0]])
        self.tnsf = imgDatasets[self.dsNames[0]].GetGeoTransform()
        self.proj = imgDatasets[self.dsNames[0]].GetProjection()

        # Retrieve phase information from IFG data
        phsMaps = []
        for dsName in imgDatasets.keys():
            # Retrieve phase map
            phs = imgDatasets[dsName].GetRasterBand(2).ReadAsArray()

            # Replace nan values
            phs[np.isnan(phs) == 1] = 0

            # Append to list
            phsMaps.append(phs)

        # Parse geometry data
        incMaps = []
        azMaps = []
        for dsName in geomDatasets.keys():
            # Retrieve image arrays
            incMap, azMap = parse_isce_los(geomDatasets[dsName])

            # Replace nan values
            incMap[np.isnan(incMap) == 1] = 0
            azMap[np.isnan(azMap) == 1] = 0

            # Append to list
            incMaps.append(incMap)
            azMaps.append(azMap)

        # Convert to numpy arrays
        self.phs = np.array(phsMaps)
        self.inc = np.array(incMaps)
        self.az = np.array(azMaps)


    ## Masking
    def __create_mask__(self, maskArgs):
        '''
        Create a common mask for the data set. Assume only the background value.
        '''
        # Pass all values initially
        self.mask = np.ones(self.phs[0,:,:].shape)

        # Loop through maps
        for i in range(self.nObs):
            # Mask by phase
            self.mask[create_mask(self.phs[i,:,:], maskArgs) == 0] = 0

            # Mask by incidence
            self.mask[create_mask(self.inc[i,:,:], maskArgs) == 0] = 0

            # Mask by azimuth
            self.mask[create_mask(self.az[i,:,:], maskArgs) == 0] = 0


    ## Vector conversion
    def __geom_to_vectors__(self):
        '''
        Convert geometries to look vectors.
        Wrapper for aria_geom_to_vector and isce_geom_to_vector.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Converting radar geometry angles to vectors')
        print(self.inc.shape)

        # Setup
        self.Px = np.zeros(self.inc.shape)
        self.Py = np.zeros(self.inc.shape)
        self.Pz = np.zeros(self.inc.shape)

        # Convert to vectors based on convention
        if self.convention in ['aria']:
            # ... based on ARIA convention
            for i in range(self.nObs):
                # Add layer to nObs x M x N array
                self.Px[i,:,:], self.Py[i,:,:], self.Pz[i,:,:] = \
                    aria_geom_to_vector(self.inc[i,:,:], self.az[i,:,:], verbose=self.verbose)

        elif self.convention in ['isce']:
            # ... based on ISCE convention
            for i in range(self.nObs):
                # Add layer to nObs x M x N array
                self.Px[i,:,:], self.Py[i,:,:], self.Pz[i,:,:] = \
                    isce_geom_to_vector(self.inc[i,:,:], self.az[i,:,:], verbose=self.verbose)


    ## Recomposition
    def __recompose__(self):
        '''
        Recompose the signal by inverting for displacement.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Recomposing signal from LOS observations')

        # Compute displacement/velocity components
        self.Uhat = invert_from_los(self.phs, self.Px, self.Py, self.Pz,
            mask=self.mask,
            nComponents=self.nComponents,
            maxCPUs=self.maxCPUs,
            verbose=self.verbose)

        # Reformat components into list
        self.Uhat = [self.Uhat[i,:,:] for i in range(self.nComponents)]


    ## Saving
    def save(self, outName):
        '''
        Save the vector components fields to a three-band GeoTiff using GDAL.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Saving to file')

        # Check outname
        confirm_outdir(outName)  # confirm output directory exists
        outName = confirm_outname_ext(outName, ['tif', 'tiff'])  # confirm output extension if GeoTiff

        # Save to GDAL data set
        save_gdal_dataset(outName, self.Uhat, mask=self.mask,
            proj=self.proj, tnsf=self.tnsf, verbose=self.verbose)


    ## Plotting
    def plot_inputs(self):
        '''
        Plot input data.
        '''
        cbarOrient = 'horizontal'

        # Loop through input files
        for i in range(self.nObs):
            # Spawn fresh figure and axes
            fig, [axPhs, axInc, axAz] = plt.subplots(ncols=3)

            # Plot phase
            fig, axPhs = plot_raster(self.phs[i,:,:], mask=self.mask, extent=self.extent,
                cmap='viridis', minPct=1, maxPct=99, cbarOrient=cbarOrient,
                fig=fig, ax=axPhs)
            axPhs.set_title('Displacement/Velocity')

            fig, axInc = plot_raster(self.inc[i,:,:], mask=self.mask, extent=self.extent,
                cmap='viridis', cbarOrient=cbarOrient,
                fig=fig, ax=axInc)
            axInc.set_title('Incidence')

            fig, axAz = plot_raster(self.az[i,:,:], mask=self.mask, extent=self.extent,
                cmap='viridis', cbarOrient=cbarOrient,
                fig=fig, ax=axAz)
            axAz.set_title('Azimuth')

            # Format figure
            fig.suptitle(self.dsNames[i])
            fig.tight_layout()


    def plot_outputs(self):
        '''
        Plot output solution.
        '''
        cbarOrient = 'horizontal'

        # Labels
        if self.nComponents == 2: titles = ['east', 'vertical']
        if self.nComponents == 3: titles = ['east', 'north', 'vertical']

        # Plot displacement/velocity componenent maps
        raster_multiplot(self.Uhat, ncols=self.nComponents,
            mask=self.mask, extent=self.extent,
            cmap='viridis', cbarOrient=cbarOrient,
            minPct=1, maxPct=99,
            titles=titles)



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Recompose signal
    recomp = LOSrecomposition(convention=inps.convention,
        imgFiles=inps.imgFiles,
        incFiles=inps.incFiles, azFiles=inps.azFiles, geomFiles=inps.geomFiles,
        maskArgs=inps.maskArgs,
        nComponents=inps.nComponents,
        maxCPUs=inps.maxCPUs,
        verbose=inps.verbose)

    # Plot inputs if requested
    if inps.plotInputs == True: recomp.plot_inputs()

    # Save to file
    recomp.save(inps.outName)

    # Plot outputs if requested
    if inps.plot == True: recomp.plot_outputs()


    plt.show()