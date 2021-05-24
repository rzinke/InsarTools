#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Read a MintPy timeseries.h5 file. Plot the points and fit a curve to them.

FUTURE IMPROVEMENTS
    * Pixel width greater than 1

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
from IOsupport import confirm_outdir, load_gdal_dataset
from GeoFormatting import transform_to_extent, lola_to_xy
from Masking import create_mask
from Fitting import dates_to_datetimes, time_since_reference, fit_linear, fit_periodic


### PARSER ---
Description = '''Read a MintPy timeseries.h5 file. Plot the points and fit a curve to them.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='avirisName', type=str,
        help='AVIRIS _corr file.')
    InputArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')


    QueryArgs = parser.add_argument_group('QUERY ARGUMENTS')
    QueryArgs.add_argument('-q','--queryfile', dest='queryFile', type=str, default=None,
        help='Query points with points in geographic coordinates, one query point per line (e.g., -118.325,35.456).')

    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name ([Out]).')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### TS ANALYSIS CLASS ---
class AVIRISanalysis:
    def __init__(self, avirisName, outName, maskArgs=None, verbose=False):
        '''
        Evaluate an AVIRIS spectral data set.
        '''
        # Parameters
        self.k = 0  # start query counter
        self.verbose = verbose

        # Load AVIRIS data set using GDAL
        self.__load_data__(avirisName)

        # Create mask
        self.mask = create_mask(self.bands[4,:,:], maskArgs, verbose=self.verbose)

        # Keep table of queried data points
        self.queriedPoints = {}

        # Format outName
        self.__format_outname__(outName)


    def __load_data__(self, avirisName):
        '''
        Load the AVIRIS data as GDAL data set.
        '''
        # Load data set
        DS = load_gdal_dataset(avirisName, verbose=self.verbose)

        # Loop though bands
        self.nBands = DS.RasterCount  # number of bands
        self.wavelengths = np.linspace(400, 2500, self.nBands)
        self.bands = []  # empty list of bands

        if self.verbose == True: print('Loading {:d} bands'.format(self.nBands))

        self.bands = [DS.GetRasterBand(i).ReadAsArray() for i in range(1, self.nBands+1)]

        self.bands = np.array(self.bands)

        # Data set sizes
        self.M, self.N = DS.RasterYSize, DS.RasterXSize

        # Geographic transform
        self.tnsf = DS.GetGeoTransform()

        # Geographic extent
        self.extent = transform_to_extent(self.tnsf, self.M, self.N, verbose=self.verbose)


    def __format_outname__(self, outName):
        '''
        Format the output name.
        '''
        confirm_outdir(outName)

        self.outName = outName


    ## Timeseries analysis
    def __get_spectrum__(self, lon, lat):
        '''
        Retrieve spectral values for the given pixel value.
        '''
        # Convert to pixel values
        px, py = lola_to_xy(self.tnsf, lon, lat, verbose=self.verbose)

        # Check that value is not within mask
        if self.mask[py, px] == 0:
            print('Point ({:.6f}, {:.6f}) is masked. Select another.'.format(lon, lat))
            spectrum = None

        else:
            # Update counter
            self.k += 1

            # Update queried points list
            self.queriedPoints[self.k] = {'lat': lat, 'lon': lon, 'px': px, 'py': py}

            # Retrieve timeseries
            spectrum = self.__retrieve_spectrum__(px, py)

        return spectrum


    def __retrieve_spectrum__(self, px, py):
        '''
        Retrieve spectrum from data cube.
        '''
        spectrum = self.bands[:, py, px]

        return spectrum


    ## Query from file
    def query_from_file(self, queryFile):
        '''
        Load lat, lon points from a text file to be used as query points.
        '''
        # Load lats, lons from file
        lons, lats = self.__load_query_file__(queryFile)
        nPts = len(lats)

        # Loop through each coordinate to retrieve the timeseries there
        for i in range(nPts):
            self.__get_spectrum__(lons[i], lats[i])

        # Save timeseries files
        self.save_spectra()


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
        lons = []; lats = []

        for line in lines:
            # Split text string
            coords = line.split(',')

            # Format lat/lon
            coords = [float(coord.strip(',').strip(' ').strip('\n')) for coord in coords]

            # Append lists
            lons.append(coords[0])
            lats.append(coords[1])

        return lons, lats


    ## Plotting
    def plot(self):
        '''
        Plot spectra.
        '''
        if self.verbose == True: print('Plotting spectra')

        # Spawn figure
        fig, ax = plt.subplots()

        # Loop through each spectrum in the queries list
        for ptNb in self.queriedPoints.keys():
            # Recall pixel location
            queriedPoint = self.queriedPoints[ptNb]  # recall each point
            px = queriedPoint['px']; py = queriedPoint['py']
            lon = queriedPoint['lon']; lat = queriedPoint['lat']

            # Recall spectrum
            spectrum = self.__retrieve_spectrum__(px, py)

            # Plot spectrum
            ax.semilogy(self.wavelengths, spectrum, label='{:d}'.format(ptNb))

        # Format plot
        ax.legend()
        ax.set_xlabel('wavelength (nm)')
        ax.set_ylabel('reflectance')


    ## Saving
    def __save_clicked_spectra__(self, event):
        '''
        Trigger saving all current selected profiles.
        '''
        # Use save_profiles function
        self.save_spectra()


    def save_spectra(self):
        '''
        Recall all queried profiles and save to text files.
        '''
        if self.verbose == True: print('Saving profiles')

        # Setup
        legend = '# Lat, Lon, px, py\n'
        metadata = '{:.6f}, {:.6f}, {:d}, {:d}\n'
        dataFmt = '{:.6f} {:.6f}\n'

        # Loop through each profile in the queries list
        for ptNb in self.queriedPoints.keys():
            # Recall pixel location
            queriedPoint = self.queriedPoints[ptNb]  # recall each point
            px = queriedPoint['px']; py = queriedPoint['py']
            lon = queriedPoint['lon']; lat = queriedPoint['lat']

            if self.verbose == True: print('Saving point ({:.6f}, {:.6f})'.format(lon, lat))

            # Format file name
            saveName = '{:s}_{:d}.txt'.format(self.outName, ptNb)

            # Recall the spectrum from that point
            spectrum = self.__retrieve_spectrum__(px, py)

            # Open file
            with open(saveName, 'w') as outFile:
                # Write legend and metadata to header
                outFile.write(legend)
                outFile.write(metadata.format(lat, lon, px, py))
                outFile.write('wavelength intensity\n')

                # Write timeseries data to file
                [outFile.write(dataFmt.format(self.wavelengths[i], spectrum[i])) for i in range(self.nBands)]

            # Report if requested
            if self.verbose == True: print('Spectrum saved to: {:s}'.format(saveName))

        # Report when finished
        if self.verbose == True: print('Saving complete')



### MAIN ---
if __name__ == '__main__':
    ## Setup
    # Gather arguments
    inps = cmdParser()


    ## 1D timeseries
    AV = AVIRISanalysis(avirisName=inps.avirisName, outName=inps.outName, maskArgs=inps.maskArgs,
        verbose=inps.verbose)


    ## Query points
    if inps.queryFile is not None:
        AV.query_from_file(inps.queryFile)


    ## Plot
    if inps.plot == True:
        AV.plot()
        plt.show()