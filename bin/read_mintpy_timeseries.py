#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Read a MintPy timeseries.h5 file. Plot the points and fit a curve to them.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from pandas.plotting import register_matplotlib_converters
import h5py
from IOsupport import confirm_outdir, load_mintpy_timeseries, load_mintpy_velocity, load_gdal_dataset
from GeoFormatting import get_mintpy_transform, transform_to_extent, lola_to_xy, get_mintpy_reference_point
from Masking import create_mask
from Viewing import plot_raster
from Fitting import dates_to_datetimes, time_since_reference, fit_linear, fit_periodic


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


    QueryArgs = parser.add_argument_group('QUERY ARGUMENTS')
    QueryArgs.add_argument('-q','--queryfile', dest='queryFile', type=str, default=None,
        help='Query points with points in geographic coordinates, one query point per line.')
    QueryArgs.add_argument('-i','--interactive', dest='interactive', action='store_true',
        help='Interactive mode.')
    QueryArgs.add_argument('--display-map', dest='displayMap', type=str, default=None,
        help='Map layer to display. Provide name of velocity file as .tif or .h5. ([None] uses cumulative displacement).')
    QueryArgs.add_argument('--fit-type', dest='fitType', type=str, default=None,
        help='Fit function ([None], linear, seasonal).')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### LOAD PLOTTING MAP ---
class MintPyTSanalysis:
    def __init__(self, tsName, outName, maskArgs=None, fitType=None, verbose=False):
        '''
        Evaluate a displacement timeseries.

        INPUTS
            dates is a E-length list of dates
            disps is a (E x M x N) array of maps, where E is the number of
             displacement epochs

        INHERITANCES
        '''
        # Parameters
        self.k = 0  # start query counter
        self.verbose = verbose

        # Format fit type
        self.__format_fit_type__(fitType)

        # Load MintPy data set
        self.__load_data__(tsName)

        # Create mask
        self.mask = create_mask(self.disps[-1,:,:], maskArgs, verbose=self.verbose)

        # Retrieve spatial extent
        self.__parse_spatial_info__(tsName)

        # Keep table of queried data points
        self.queriedPoints = {}

        # Format outName
        self.__format_outname__(outName)


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
        self.E, self.M, self.N = self.disps.shape

        # Geographic transform
        self.tnsf = get_mintpy_transform(tsName, verbose=inps.verbose)

        # Geographic extent
        self.extent = transform_to_extent(self.tnsf, self.M, self.N, verbose=self.verbose)

        # Reference point
        self.refLon, self.refLat = get_mintpy_reference_point(tsName, verbose=self.verbose)


    def __format_outname__(self, outName):
        '''
        Format the output name.
        '''
        confirm_outdir(outName)

        self.outName = outName


    ## Timeseries analysis
    def __get_displacements__(self, lon, lat):
        '''
        Retrieve epoch displacements for the given pixel value.
        '''
        # Convert to pixel values
        px, py = lola_to_xy(self.tnsf, lon, lat, verbose=self.verbose)

        # Check that value is not within mask
        if self.mask[py, px] == 0:
            print('Point ({:.6f}, {:.6f}) is masked. Select another.'.format(lat, lon))
            ts = None

        else:
            # Update counter
            self.k += 1

            # Update queried points list
            self.queriedPoints[self.k] = {'lat': lat, 'lon': lon, 'px': px, 'py': py}

            # Retrieve timeseries
            ts = self.__retrieve_ts__(px, py)

        return ts


    def __retrieve_ts__(self, px, py):
        '''
        Retrieve timeseries from data cube.
        '''
        ts = self.disps[:, py, px]

        return ts


    ## Query from file
    def query_from_file(self, queryFile):
        '''
        Load lat, lon points from a text file to be used as query points.
        '''
        # Load lats, lons from file
        lats, lons = self.__load_query_file__(queryFile)
        nPts = len(lats)

        # Loop through each coordinate to retrieve the timeseries there
        for i in range(nPts):
            self.__get_displacements__(lons[i], lats[i])

        # Save timeseries files
        self.save_profiles()


    def __load_query_file__(self, queryFile):
        '''
        Load points from the specified query file.
        '''
        # Open and load file contents
        with open(queryFile, 'r') as qFile:
            # Get lines
            lines = qFile.readlines()
            lines = [line for line in lines if line[0]!='#']  # ignore header

        # Convert to lats, lons
        lats = []; lons = []

        for line in lines:
            # Split text string
            coords = line.split()

            # Format lat/lon
            coords = [float(coord.strip(',').strip('\n')) for coord in coords]

            # Append lists
            lats.append(coords[0])
            lons.append(coords[1])

        return lats, lons


    ## Interactive plotting
    def interactive(self, basemapName=None):
        '''
        Interactive plot for timeseries analysis.
        '''
        # Determine base map
        if basemapName is not None:
            self.basemap = self.__load_base_map__(basemapName)
        else:
            self.basemap = self.disps[-1,:,:]  # final displacement

        # Establish timeseries plot
        self.tsFig = plt.figure(figsize=(10, 6))

        # Timeseries axis
        self.axTS = self.tsFig.add_subplot(position=(0.1, 0.3, 0.8, 0.6))

        # Clear profiles button
        self.axClear = self.tsFig.add_subplot(position=(0.15, 0.05, 0.1, 0.05))
        self.clearButton = Button(self.axClear, 'Clear profile', hovercolor='1')
        self.clearButton.on_clicked(self.__clear_figs__)

        # Save profile button
        self.axSave = self.tsFig.add_subplot(position=(0.45, 0.05, 0.1, 0.05))
        self.saveButton = Button(self.axSave, 'Save profiles', hovercolor='1')
        self.saveButton.on_clicked(self.__save_clicked_profiles__)

        # Establish map plot
        self.mapFig, self.axMap = plt.subplots(figsize=(5, 8))
        self.__plot_basemap__()

        # Interact with map
        self.mapFig.canvas.mpl_connect('button_press_event', self.__plot_point_ts__)


    def __load_base_map__(basemapName=None):
        '''
        Load the specified map file. Automatically detect the file type based on 
         the file extension.
        '''
        # Load the data set according to the file type
        if self.verbose == True: print('Loading map layer for plotting')

        # Idenitfy file type
        ext = basemapName.split('.')[-1]  # retrieve extension
        assert ext in ['tif', 'h5'], 'Cannot identify display map extension {:s}'.format(ext)

        # Load file based on extension
        if ext == 'tif':
            # Load map as geotiff
            DS = load_gdal_dataset(basemapName, verbose=self.verbose)
            self.plotLayer = DS.GetRasterBand(1).ReadAsArray()
        elif ext == 'h5':
            # Load map as MintPy HDF5 data set - assume that layer is velocity
            self.plotLayer = load_mintpy_velocity(basemapName, verbose=self.verbose)


    def __plot_basemap__(self):
        '''
        Plot base map.
        '''
        # Plot base map
        self.mapFig, self.axMap = plot_raster(self.basemap,
            mask=self.mask, extent=self.extent,
            minPct=1, maxPct=99, cbarOrient='vertical',
            fig=self.mapFig, ax=self.axMap)

        # Plot reference point
        self.axMap.plot(self.refLon, self.refLat, 'ks')


    def __plot_point_ts__(self, event):
        '''
        Plot the displacement history of a single point.
        '''
        # Record click values
        lon = event.xdata
        lat = event.ydata

        # Retrieve displacement values
        ts = self.__get_displacements__(lon, lat)

        # Check that value is not within mask
        if ts is not None:
            # Random color
            self.color = np.random.rand(3)  # random color

            # Plot query point on map
            self.axMap.scatter(lon, lat, facecolor=self.color, edgecolor='k')

            # Plot timeseries
            self.axTS.scatter(self.datetimes, ts, s=25, c=[self.color]*len(ts))

            # Fit timeseries
            self.__fit_timeseries__(ts)

            # Render images
            self.mapFig.canvas.draw()
            self.tsFig.canvas.draw()


    def __fit_timeseries__(self, ts):
        '''
        Fit a timeseries according to the prescribed fit type.
        '''
        # Format label
        label='Profile {:d}'.format(self.k)

        # Fit
        if self.fitType is not None:
            # Fit by fit type
            if self.fitType == 'linear':
                # Fit linear
                fit, B = fit_linear(self.times, ts, verbose=self.verbose)
                self.axTS.plot(self.datetimes, fit, color=self.color, label=label)

            elif self.fitType == 'seasonal':
                # Fit seasonal
                fit, B = fit_periodic(self.times, ts, verbose=self.verbose)
                self.axTS.plot(self.datetimes, fit, color=self.color, label=label)

        # No fit
        else:
            self.axTS.plot(self.datetimes[0], ts[0], color=self.color, label=label)

        # Update legend
        self.axTS.legend()


    def __clear_figs__(self, event):
        '''
        Clear lists and figures.
        '''
        # Clear queried points list
        self.queriedPoints = {}

        # Clear map figure
        self.axMap.images[-1].colorbar.remove()
        self.axMap.cla()

        # Replot map
        self.__plot_basemap__()

        # Clear timeseries fig
        self.axTS.cla()

        # Render images
        self.mapFig.canvas.draw()
        self.tsFig.canvas.draw()

        if self.verbose == True: print('Figure and cache cleared')


    ## Saving
    def __save_clicked_profiles__(self, event):
        '''
        Trigger saving all current selected profiles.
        '''
        # Use save_profiles function
        self.save_profiles()


    def save_profiles(self):
        '''
        Recall all queried profiles and save to text files.
        '''
        if self.verbose == True: print('Saving profiles')

        # Setup
        legend = '# Lat, Lon, px, py / Date (YYYYMMDD), displacement (m)\n'
        metadata = '{:.6f}, {:.6f}, {:d}, {:d}\n'
        dataFmt = '{:s}, {:.6f}\n'

        # Loop through each profile in the queries list
        for ptNb in self.queriedPoints.keys():
            # Recall pixel location
            queriedPoint = self.queriedPoints[ptNb]  # recall each point
            px = queriedPoint['px']; py = queriedPoint['py']
            lon = queriedPoint['lon']; lat = queriedPoint['lat']

            if self.verbose == True: print('Saving point ({:.6f}, {:.6f})'.format(lat, lon))

            # Format file name
            saveName = '{:s}_{:d}.txt'.format(self.outName, ptNb)

            # Recall the displacements from that point
            ts = self.__retrieve_ts__(px, py)

            # Open file
            with open(saveName, 'w') as outFile:
                # Write legend and metadata to header
                outFile.write(legend)
                outFile.write(metadata.format(lat, lon, px, py))

                # Write timeseries data to file
                [outFile.write(dataFmt.format(self.dates[i], ts[i])) for i in range(len(ts))]

            # Report if requested
            if self.verbose == True: print('Profile saved to: {:s}'.format(saveName))

        # Report when finished
        if self.verbose == True: print('Saving complete')



### MAIN ---
if __name__ == '__main__':
    ## Setup
    # Gather arguments
    inps = cmdParser()


    ## 1D timeseries
    TS = MintPyTSanalysis(tsName=inps.tsName, outName=inps.outName, maskArgs=inps.maskArgs, fitType=inps.fitType,
        verbose=inps.verbose)


    ## Query points
    if inps.queryFile is not None:
        TS.query_from_file(inps.queryFile)


    ## Interactive mode
    # Plot map
    if inps.interactive == True:
        TS.interactive(inps.displayMap)
        plt.show()