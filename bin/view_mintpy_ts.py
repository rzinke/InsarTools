#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Read a MintPy timeseries.h5 file. Plot the points and fit a curve to them.

FUTURE IMPROVEMENTS

TESTING STATUS
In development.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters
import h5py
from IOsupport import load_mintpy_timeseries
from GeoFormatting import get_mintpy_transform, transform_to_extent, get_mintpy_reference_point
from Fitting import dates_to_datetimes, time_since_reference
from Masking import create_mask
from Viewing import plot_raster


### PARSER ---
Description = '''Read a MintPy timeseries.h5 file. Plot the points and fit a curve to them.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='tsName', type=str,
        help='MintPy timeseries.h5 file.')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')


    DisplayArgs = parser.add_argument_group('DISPLAY PARAMS')
    DisplayArgs.add_argument('-i','--start-index', dest='startNdx', type=int, default=0,
        help='First index to plot.')
    DisplayArgs.add_argument('-d','--start-date', dest='startDate', type=str, default=None,
        help='First date to plot.')
    DisplayArgs.add_argument('-c','--cmap', dest='cmap', type=str, default='viridis',
        help='Colormap ([viridis]).')
    DisplayArgs.add_argument('-co','--colorbar-orientation', dest='cbarOrient', type=str, default='horizontal',
        help='Colorbar orientation ([horizontal], vertical).')
    DisplayArgs.add_argument('-minPct','--min-percent', dest='minPct', type=float, default=None,
        help='Minimum percent clip value ([None]).')
    DisplayArgs.add_argument('-maxPct','--max-percent', dest='maxPct', type=float, default=None,
        help='Maximum percent clip value ([None]).')
    DisplayArgs.add_argument('-vmin','--min-value', dest='vmin', type=float, default=None,
        help='Minimum clip value ([None]).')
    DisplayArgs.add_argument('-vmax','--max-value', dest='vmax', type=float, default=None,
        help='Maximum clip value ([None]).')
    DisplayArgs.add_argument('--steady-colors', dest='steadyColors', action='store_true',
        help='Use the min/max value of the final displacement map, clipped at 1 percent.')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### TS ANALYSIS CLASS ---
class MintPyTSview:
    def __init__(self, tsName, maskArgs=None, verbose=False):
        '''
        Evaluate a displacement timeseries.

        INPUTS
            dates is a E-length list of dates
            disps is a (E x M x N) array of maps, where E is the number of
             displacement epochs
        '''
        # Parameters
        self.verbose = verbose

        # Load MintPy data set
        self.__load_data__(tsName)

        # Create mask
        self.mask = create_mask(self.disps[-1,:,:], maskArgs, verbose=self.verbose)

        # Retrieve spatial extent
        self.__parse_spatial_info__(tsName)

        # Show dates
        if self.verbose == True:
            print('Available epochs:')
            [print('{:d} {:s}'.format(i, date)) for i, date in enumerate(self.dates)]


    def __load_data__(self, tsName):
        '''
        Load the MintPy (HDF5) data set, parse the dates and displacements.
        '''
        # Open and close HDF5 data set
        self.dates, self.disps = load_mintpy_timeseries(tsName, verbose=self.verbose)

        # Convert dates to datetimes
        self.datetimes = dates_to_datetimes(self.dates, verbose=self.verbose)

        # Calculate time since beginning of series
        self.times = np.array(time_since_reference(self.datetimes, verbose=self.verbose))


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


    def plot(self, startNdx=0, startDate=None, cmap='viridis', cbarOrient='vertical',
        vmin=None, vmax=None, minPct=1, maxPct=99,
        steadyColors=False):
        '''
        Initial plot of the timeseries.
        '''
        # Parameters
        self.cmap = cmap
        self.cbarOrient = cbarOrient
        self.vmin = vmin
        self.vmax = vmax
        self.minPct = minPct
        self.maxPct = maxPct

        # Define min/max colors or "steady colors" is specified
        if steadyColors == True:
            self.__compute_minmax_colors__()

        # Spawn figure
        self.fig, self.ax = plt.subplots()

        # Determine index to plot
        self.ndx = self.__find_start_date__(startDate)

        # Plot initial displacement to start
        self.__plot_epoch__()

        # Connect to keyboard
        self.fig.canvas.mpl_connect('key_press_event', self.__switch_maps__)


    def __compute_minmax_colors__(self):
        '''
        Use the min/max value of the final displacement map, clipped at
         1 percent.
        '''
        # Determine clip values
        Dfinal = np.ma.array(self.disps[-1,:,:], mask=(self.mask==0))
        Dfinal = Dfinal.compressed().flatten()
        self.vmin, self.vmax = np.percentile(Dfinal, [1, 99])

        # Report if requested
        if self.verbose == True: print('Clipping colors to: {:f}, {:f}'.format(self.vmin, self.vmax))


    def __find_start_date__(self, startDate):
        '''
        Find the index of the given start date in YYYYMMDD form.
        '''
        if startDate is not None:
            ndx = self.dates.index(startDate)

            if self.verbose == True:
                print('Displaying {:d} {:s}'.format(ndx, startDate))
        else:
            ndx = 0

        return ndx


    def __plot_epoch__(self):
        '''
        Plot an epoch by the index provided.
        '''
        # Plot raster
        self.fig, self.ax = plot_raster(self.disps[self.ndx,:,:],
            mask=self.mask, extent=self.extent,
            cmap=self.cmap, cbarOrient=self.cbarOrient,
            vmin=self.vmin, vmax=self.vmax,
            minPct=self.minPct, maxPct=self.maxPct,
            fig=self.fig, ax=self.ax)

        # Format plot
        self.ax.set_title(self.dates[self.ndx])


    def __switch_maps__(self, event):
        '''
        Switch maps foward or backward depending on the key press.
        '''
        # Clear axis
        self.ax.set_title(None)  # clear title
        self.ax.images[-1].colorbar.remove()  # clear colorbar
        self.ax.cla()  # clear map

        # Register left or right key press
        if event.key == 'right':
            # Move forward in time
            self.__update_index__(1)  # update figure index
            self.__plot_epoch__()

        elif event.key == 'left':
            # Move backward in time
            self.__update_index__(-1)  # update figure index
            self.__plot_epoch__()

        else:
            print('Use <-left or right-> arrow keys to move through data set.')

        # Render plot
        self.fig.canvas.draw()


    def __update_index__(self, increment):
        '''
        Update the epoch index, if that is valid.
        '''
        # Check if update is valid
        if (self.ndx + increment < self.nEpochs) and (self.ndx + increment >= 0):
            # Update index
            self.ndx += increment



### MAIN ---
if __name__ == '__main__':
    ## Setup
    # Gather arguments
    inps = cmdParser()


    ## Plot timeseries
    # Load data set
    TS = MintPyTSview(tsName=inps.tsName, maskArgs=inps.maskArgs, verbose=inps.verbose)

    # Plot
    TS.plot(startNdx=inps.startNdx, startDate=inps.startDate,
        cmap=inps.cmap, cbarOrient=inps.cbarOrient,
        vmin=inps.vmin, vmax=inps.vmax,
        minPct=inps.minPct, maxPct=inps.maxPct,
        steadyColors=inps.steadyColors)


    print('__In development__')


    plt.show()