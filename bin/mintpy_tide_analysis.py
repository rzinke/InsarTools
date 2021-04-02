#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Fit a ramp to the epoch displacements in the MintPy timeseries, and the corresponding solid Earth tide measurements.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
import h5py
from IOsupport import confirm_outdir, confirm_outname_ext, load_mintpy_timeseries
from GeoFormatting import get_mintpy_transform, transform_to_extent, lola_to_xy, get_mintpy_reference_point, get_mintpy_reference_pixel
from Masking import create_mask
from Viewing import plot_raster
from Fitting import dates_to_datetimes, time_since_reference, fit_surface


### PARSER ---
Description = '''Fit a ramp to the epoch displacements in the MintPy timeseries, and the corresponding solid Earth tide measurements.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument('-t','--timeseries', dest='tsName', type=str,
        help='MintPy timeseries.h5 file.')
    InputArgs.add_argument('-s','--solidtide', dest='setName', type=str,
        help='MintPy SET.h5 file.')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')
    InputArgs.add_argument('-d','--decimation-factor', dest='decimationFactor', type=float, default=0,
        help='Decimation factor for fitting.')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot series.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### TS ANALYSIS CLASS ---
class MintPyTideAnalysis:
    def __init__(self, tsName, setName, outName, maskArgs=None, decimationFactor=0, verbose=False):
        '''
        Fit a ramp to the epochs in the MintPy displacement timeseries and the
         corresponding solid earth tide measurements.
        '''
        # Parameters
        self.verbose = verbose

        # Load MintPy data set
        self.__load_data__(tsName, setName)

        # Create mask
        self.mask = create_mask(self.disps[-1,:,:], maskArgs, verbose=self.verbose)

        # Retrieve spatial extent
        self.__parse_spatial_info__(tsName)

        # Fit plane to displacements
        self.__fit_displacements__(decimationFactor)

        # Format outName
        self.__format_outname__(outName)

        # Save to file
        self.__save__()


    def __load_data__(self, tsName, setName):
        '''
        Load the MintPy (HDF5) data set, parse the dates and displacements.
        '''
        # Open and close HDF5 data set
        self.dates, self.disps = load_mintpy_timeseries(tsName, verbose=self.verbose)

        # Convert dates to datetimes
        self.datetimes = dates_to_datetimes(self.dates, verbose=self.verbose)

        # Calculate time since beginning of series
        self.times = np.array(time_since_reference(self.datetimes, verbose=self.verbose))

        # Load solid Earth tide data
        _, self.tides = load_mintpy_timeseries(setName, verbose=self.verbose)


    def __parse_spatial_info__(self, tsName):
        '''
        Retrieve the geographic information
        '''
        # Data set sizes
        self.nEpochs, self.M, self.N = self.disps.shape

        # Geographic transform
        self.tnsf = get_mintpy_transform(tsName, verbose=inps.verbose)

        # Geographic extent
        self.extent = transform_to_extent(self.tnsf, self.M, self.N, verbose=self.verbose)

        # Reference point
        self.refLon, self.refLat = get_mintpy_reference_point(tsName, verbose=self.verbose)
        self.refX, self.refY = get_mintpy_reference_pixel(tsName, verbose=self.verbose)


    def __fit_displacements__(self, decimationFactor):
        '''
        Fit a plane to the displacement timeseries.
        '''
        # Setup
        self.dispFits = []
        self.tideFits = []

        # Loop through displacements
        for i in range(3): #range(self.nEpochs):
            if self.verbose == True: print('Epoch: {:s}'.format(self.dates[i]))
            # Fit surface to displacements
            dispSurface, B = fit_surface(img=self.disps[i], mask=self.mask,
                decimation=decimationFactor,
                verbose=self.verbose)

            fig, axes = plt.subplots(figsize=(8,8), nrows=2, ncols=3)

            disps = np.ma.array(self.disps[i], mask=(self.mask==0))
            vmin, vmax = disps.compressed().min(), disps.compressed().max()
            axes[0][0].imshow(disps)
            axes[0][1].imshow(dispSurface)
            axes[0][2].imshow(disps-dispSurface, vmin=vmin, vmax=vmax)

            # Append to list
            self.dispFits.append(B)


            # Fit surface to solid Earth tides
            tideSurface, B = fit_surface(img=self.tides[i], mask=self.mask,
                decimation=decimationFactor,
                verbose=self.verbose)

            tide = np.ma.array(self.tides[i], mask=(self.mask==0))
            vmin, vmax = tide.compressed().min(), tide.compressed().max()
            axes[1][0].imshow(tide)
            axes[1][1].imshow(tideSurface)
            axes[1][2].imshow(tide-tideSurface, vmin=vmin, vmax=vmax)

            fig.savefig('TideComparisonsRefPt/{:s}.png'.format(self.dates[i]))

            # Append to list
            self.tideFits.append(B)

        # plt.show()
        exit()


    def __format_outname__(self, outName):
        '''
        Format the output name.
        '''
        outName = confirm_outname_ext(outName, 'txt')
        confirm_outdir(outName)

        self.outName = outName


    def __save__(self):
        '''
        Save results to file.
        '''
        with open(self.outName, 'w') as outFile:
            outFile.write('# date, fit coefficients\n')
            for i in range(self.nEpochs):
                outFile.write('{}, {}, {}\n'.format(self.dates[i], self.dispFits[i], self.tideFits[i]))


    def plot(self):
        '''
        Plot results.
        '''
        # Convert displacement fits to numpy array
        xDisps = [dispFit[1] for dispFit in self.dispFits]
        xTides = [tideFit[1] for tideFit in self.tideFits]

        # Spawn figure
        fig, ax = plt.subplots()

        ax.plot(self.times, xDisps, '-ro')
        ax.plot(self.times, xTides, '-bo')



### MAIN ---
if __name__ == '__main__':
    ## Setup
    # Gather arguments
    inps = cmdParser()


    ## Timeseries analysis
    TS = MintPyTideAnalysis(tsName=inps.tsName, setName=inps.setName, outName=inps.outName,
        maskArgs=inps.maskArgs,
        decimationFactor=inps.decimationFactor,
        verbose=inps.verbose)


    if inps.plot == True:
        TS.plot()

        plt.show()