#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Convert timeseries displacements to velocities.

FUTURE IMPROVEMENTS
    * Accept list of fnames
    * Accept list of dates in command line

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from pandas.plotting import register_matplotlib_converters
import multiprocessing as mp
from IOsupport import confirm_outdir, confirm_outname_ext, append_fname, load_gdal_dataset, save_gdal_dataset
from GeoFormatting import transform_to_extent, lola_to_xy
from Masking import create_mask
from Viewing import plot_raster
from Fitting import dates_to_datetimes, time_since_reference, fit_linear, fit_periodic


### PARSER ---
Description = '''Convert timeseries displacements to velocities.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='tsName', type=str,
        help='Timeseries file(s).')
    InputArgs.add_argument('-d','--date-list', dest='dateList', type=str, required=True,
        help='List of dates in format YYYYMMDD (list or file with one date per line).')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')
    InputArgs.add_argument('-f','--fit-type', dest='fitType', type=str, default='linear',
        help='Fit type used to compute velocities.')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot mode.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### TS TO VELOCITY ---
class TS2velocity:
    def __init__(self, tsName, dateList, outName, maskArgs=None, verbose=False):
        '''
        Load a displacement timeseries and compute the velocity at each pixel location.

        INPUTS
            dates is a E-length list of dates
            disps is a (E x M x N) array of maps, where E is the number of
             displacement epochs
        '''
        # Parameters
        self.verbose = verbose

        # Load MintPy data set
        self.__load_images__(tsName)

        # Load dates
        self.__load_dates__(dateList)

        # Create mask
        self.mask = create_mask(self.disps[-1,:,:], maskArgs, verbose=self.verbose)

        # Format outName
        self.__format_outname__(outName)


    def __load_images__(self, tsName):
        '''
        Load images and parse the spatial data.
        '''
        # Load GDAL data sets
        DS = load_gdal_dataset(tsName, verbose=self.verbose)

        # Number of epochs
        self.nEpochs = DS.RasterCount

        # Format displacement maps into array
        disps = []
        for i in range(self.nEpochs):
            disps.append(DS.GetRasterBand(i+1).ReadAsArray())

        # Store as attribute
        self.disps = []
        for disp in disps:
            # Replace NaNs with 0s
            disp[np.isnan(disp)==1] = 0

            # Append to list
            self.disps.append(disp)

        self.disps = np.array(self.disps)

        # Geographic information
        self.M, self.N = DS.RasterYSize, DS.RasterXSize
        self.tnsf = DS.GetGeoTransform()
        self.proj = DS.GetProjection()

        # Report if requested
        if self.verbose == True:
            print('Number of epochs detected: {:d}'.format(self.nEpochs))
            print('Image size: {:d} x {:d}'.format(self.M, self.N))

        # Geographic extent
        self.extent = transform_to_extent(self.tnsf, self.M, self.N, verbose=self.verbose)


    def __load_dates__(self, dateList):
        '''
        Load and format list of dates.
        '''
        # Load dates from date list
        with open(dateList, 'r') as dateFile:
            dates = dateFile.readlines()

        # Format dates
        self.dates = [date.strip(',').strip('\n') for date in dates]

        # Check that number of dates is the same as number of epochs
        assert len(self.dates) == self.nEpochs, 'Number of dates must equal number of images.'

        # Convert dates to datetimes
        self.datetimes = dates_to_datetimes(self.dates, verbose=self.verbose)

        # Calculate time since beginning of series
        self.times = np.array(time_since_reference(self.datetimes, verbose=self.verbose))


    def __format_outname__(self, outName):
        '''
        Format the output name.
        '''
        confirm_outdir(outName)

        self.outName = outName


    ## Velocity computation
    def compute_velocity(self, fitType='linear'):
        '''
        Compute the secular velocity using the given fit type.
        '''
        if self.verbose == True: print('Computing velocity with fit type: {:s}'.format(fitType))

        # Check fit type and determine function
        self.__format_fit_type__(fitType)

        # Format inversion matrix
        self.__build_design_matrix__()
        self.__invert_design_matrix__()

        # Loop through each pixel
        self.velocityMap = np.zeros((self.M, self.N))
        self.residualMap = np.zeros((self.M, self.N))

        for i in range(self.M):
            for j in range(self.N):
                # Recall displacements
                disps = self.disps[:,i,j]

                # Solve for component scaling factors
                B = self.__find_velocity_components__(disps)

                # Store linear velocity
                self.velocityMap[i,j] = B[1]  # second term is always linear component

                # Compute residuals
                residuals = disps - self.G.dot(B)

                # Store residuals
                self.residualMap[i,j] = np.sqrt(np.mean(residuals**2))


    def __format_fit_type__(self, fitType):
        '''
        Check and format the type of fit to apply, if applicable.
        '''
        if fitType is not None:
            # Ensure lower case
            fitType = fitType.lower()

            # Check that fit type is acceptable
            assert fitType in ['linear', 'seasonal'], \
                'Fit type {:s} is not valid, use \'linear\' or \'seasonal\''.format(fitType)

        self.fitType = fitType


    def __build_design_matrix__(self):
        '''
        Create the design matrix based on the fit type.
        First component is always offset; second is always linear velocity.
        '''
        # For linear fit
        if self.fitType == 'linear':
            # Formulate design matrix
            self.G = np.ones((self.nEpochs, 2))
            self.G[:,1] = self.times.flatten()

        # For seasonal fit
        elif self.fitType == 'seasonal':
            # Frequency
            freq = 1  # per year

            # Formulate design matrix
            self.G = np.ones((self.nEpochs, 4))
            self.G[:,1] = self.times.flatten()
            self.G[:,2] = np.sin(2*np.pi*freq*self.times).flatten()
            self.G[:,3] = np.cos(2*np.pi*freq*self.times).flatten()


    def __invert_design_matrix__(self):
        '''
        Invert the design matrix based on the fit type.
        Ideally only do this once.
        '''
        # Invert design matrix
        self.inverseMatrix = np.linalg.inv(np.dot(self.G.T, self.G)).dot(self.G.T)


    def __find_velocity_components__(self, displacements):
        '''
        Solve for the beta vector by taking the dot product Ginv.displacements.
        The second term in the beta vector is the linear component.

        Feed this through a loop to parallelize.
        '''
        B = self.inverseMatrix.dot(displacements)

        return B


    ## Plotting
    def plot(self):
        '''
        Plot velocity map.
        '''
        # Determine orientations
        if self.M > self.N:
            figsize = (5, 8)
            cbarOrient = 'vertical'  # colorbar orientation
        else:
            figsize = (8, 5)
            cbarOrient = 'horizontal'  # colorbar orientation

        # Velocity map
        velFig, axVel = plt.subplots(figsize=figsize)

        # Plot velocity map
        plot_raster(self.velocityMap, mask=self.mask, extent=self.extent,
            cmap='jet', cbarOrient=cbarOrient,
            minPct=1, maxPct=99, fig=velFig, ax=axVel)
        axVel.set_title('Linear velocity')

        # Residual map
        resFig, axRes = plt.subplots(figsize=figsize)

        # Plot residual map
        plot_raster(self.residualMap, mask=self.mask, extent=self.extent,
            cmap='viridis', cbarOrient=cbarOrient,
            minPct=1, maxPct=99, fig=resFig, ax=axRes)
        axRes.set_title('Residuals')


    ## Saving
    def save(self, outName):
        '''
        Save velocity map to GDAL data set.
        '''
        if self.verbose == True: print('Saving to file: {:s}'.format(outName))

        confirm_outdir(outName)  # confirm output directory exists
        velName = confirm_outname_ext(outName, ['tif', 'tiff'])  # confirm file extension

        # Save velocity map
        save_gdal_dataset(velName, self.velocityMap, mask=self.mask,
            proj=self.proj, tnsf=self.tnsf, verbose=inps.verbose)

        # Format residuals name
        resName = append_fname(velName, '_residuals')

        # Save residuals map
        save_gdal_dataset(resName, self.residualMap, mask=self.mask,
            proj=self.proj, tnsf=self.tnsf, verbose=inps.verbose)




### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather inputs
    inps = cmdParser()


    ## Load data
    TS = TS2velocity(tsName=inps.tsName, dateList=inps.dateList, outName=inps.outName, maskArgs=inps.maskArgs,
        verbose=inps.verbose)

    # Compute velocity map
    TS.compute_velocity(fitType=inps.fitType)

    # Plot velocity map
    TS.plot()

    # Save velocity map
    TS.save(inps.outName)


    plt.show()