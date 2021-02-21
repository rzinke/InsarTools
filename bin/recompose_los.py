#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Convert two or more lines of sight into movement components.

FUTURE IMPROVEMENTS

TESTING STATUS
In development. Pick up at load_data_isce and plot_inputs
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import load_gdal_datasets
from Checks import check_dataset_sizes
from RadarGeometries import parse_isce_los
from RasterResampling import match_rasters
from Viewing import plot_raster
from Masking import create_mask


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

    InputArgs.add_argument('-a','--azimuth', dest='azFiles', type=str, nargs='+', default=None,
        help='ARIA azimuth files.')
    InputArgs.add_argument('-i','--incidence', dest='incFiles', type=str, nargs='+', default=None,
        help='ARIA incidence files.')
    InputArgs.add_argument('-g','--geometry', dest='geomFiles', type=str, nargs='+', default=None,
        help='ISCE geometry files.')
    InputArgs.add_argument('-f','--ifg', dest='imgFiles', type=str, nargs='+', required=True,
        help='LOS file names.')

    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('--plot-inputs', dest='plotInputs', action='store_true',
        help='Plot input data.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true', 
        help='Plot outputs.')
    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### RECOMPOSITION CLASS ---
class LOSrecomposition:
    def __init__(self, convention, imgFiles, azFiles, incFiles, geomFiles, verbose=False):
        '''
        Use the satellite geometry files and LOS values to "recompose" the
         EW, (NS), and vertical components of motion.

        Inherits:
            IOsupport: load_gdal_datasets
            Checks: check_dataset_sizes
            RadarGeometries: parse_isce_los
            RasterResampling: match_rasters
            Masking: create_mask
            Viewing: plot_raster

        To initialize, check that inputs types are consistent.
         Then compute components of motion.
        '''
        # Parameters
        self.convention = convention.lower()
        self.verbose = verbose

        # Check that inputs are consistent
        self.__check_files__(imgFiles, azFiles, incFiles, geomFiles)

        # Load data from files
        self.__load_data__(imgFiles, azFiles, incFiles, geomFiles)

        # Mask
        self.__create_mask__()


    ## Checks
    def __check_files__(self, imgFiles, azFiles, incFiles, geomFiles):
        '''
        Check that the number of files is consistent.
        Store attribute nFiles.
        '''
        # Check number of files is consistent ...
        if self.convention in ['isce']:
            # ... for ISCE data
            assert len(imgFiles) == len(geomFiles), 'Number of image files must equal number of geometry files.'

        elif self.convention in ['aria']:
            # ... for ARIA data
            assert len(imgFiles) == len(azFiles) == (incFiles), \
                'Number of image files must equal number of azimuth and incidence files.'

        # Number of files for quick reference
        self.nFiles = len(imgFiles)


    ## Load data
    def __load_data__(self, imgFiles, azFiles, incFiles, geomFiles):
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

        # Resample data sets for consistency
        imgDatasets = match_rasters(imgDatasets, cropping='intersection', verbose=self.verbose)
        geomDatasets = match_rasters(geomDatasets, cropping='intersection', verbose=self.verbose)

        # Retrieve phase information from IFG data
        self.phsMaps = []
        for dsName in imgDatasets.keys():
            # Retrieve phase map
            phs = imgDatasets[dsName].GetRasterBand(2).ReadAsArray()

            # Replace nan values
            print(phs.shape)
            phs[phs == 'nan'] = 0

            # Append to list
            self.phsMaps.append(phs)

        # Parse geometry data
        self.incMaps = []
        self.azMaps = []
        for dsName in geomDatasets.keys():
            # Retrieve image arrays
            incMap, azMap = parse_isce_los(geomDatasets[dsName])
            self.incMaps.append(incMap)
            self.azMaps.append(azMap)


    def __load_data_aria__(self, imgFiles, azFiles, incFiles):
        '''
        Load ARIA data as GDAL data sets.
        '''
        print('ARIA data not compatible yet.')
        exit()


    ## Masking
    def __create_mask__(self):
        '''
        Create a common mask for the data set. Assume only the background value.
        '''
        # Pass all values initially
        self.mask = np.ones(self.phsMaps[0].shape)

        # Loop through phase maps
        for phs in self.phsMaps:
            mask = create_mask(phs, ['bg', 'nan'])

            # Update mask
            self.mask[mask == 0] = 0


    ## Plotting
    def plot_inputs(self):
        '''
        Plot input data.
        '''
        cbarOrient = 'horizontal'

        # Loop through input files
        for i in range(self.nFiles):
            # Spawn fresh figure and axes
            fig, [axPhs, axInc, axAz] = plt.subplots(ncols=3)

            # Plot phase
            plot_raster(self.phsMaps[i], 
                cmap='jet', minPct=1, maxPct=99, cbarOrient=cbarOrient,
                fig=fig, ax=axPhs)



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Recompose signal
    recomp = LOSrecomposition(convention=inps.convention,
        imgFiles=inps.imgFiles,
        azFiles=inps.azFiles, incFiles=inps.incFiles, geomFiles=inps.geomFiles,
        verbose=inps.verbose)

    # Plot inputs if requested
    if inps.plotInputs == True:
        recomp.plot_inputs()

    print('Do it quick. Do it well.')


    plt.show()