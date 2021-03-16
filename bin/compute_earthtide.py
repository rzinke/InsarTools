#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Use the GMT 6 program earthtide to create a solid earth tide map of the area
 described by the given georeferenced file.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import os
import matplotlib.pyplot as plt
from IOsupport import load_gdal_dataset, load_gdal_datasets
from GeoFormatting import get_raster_size, parse_transform
from Viewing import raster_multiplot


### PARSER ---
Description = '''Create a solid earth tide map of the area covered by a given map data set
using the GMT 6 program earthtide.

Format solid Earth command
e.g., "gmt earthtide -T2018-06-18T12:00:00 -GXXsolidtides_%s.grd -Ce,n,v -R80/90/30/40 -I15s"
'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument('-f','--filename', dest='dsName', type=str, required=True,
        help='Filename of file for which tides are to be calculated.')
    InputArgs.add_argument('-d','--date', dest='date', type=str, required=True,
        help='Date on which the tides are to be calculated (fmt YYYYMMDD).')
    InputArgs.add_argument('-t','--time', dest='time', type=str, required=True,
        help='Time at which the tides are to be calculated (fmt HHmmSS).')

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



### TIDE MAP ---
class create_tide_map:
    def __init__(self, dsName, date, time, outName, verbose=False):
        '''
        Create one or more tide maps using GMT 6's earthtide functionality.
        '''
        # Parameters
        self.verbose = verbose

        # Extract spatial parameters from data set
        self.__get_spatial_data__(dsName)

        # Format date and time
        self.__format_date_time__(date, time)

        # Format outname
        self.__format_outname__(outName)

        # Create tide maps
        self.__create_tide_maps__()


    def __get_spatial_data__(self, dsName):
        '''
        Retrieve the necessary spatial information (bounds and resolution) from
         the specified GDAL-compatible data set.
        '''
        if self.verbose == True: print('Retrieving data set bounds and resolution')

        # Load GDAL data set
        DS = load_gdal_dataset(dsName)
        tnsf = DS.GetGeoTransform()
        M, N = get_raster_size(DS)

        # Parse geographic info
        geo = parse_transform(tnsf, M, N, verbose=inps.verbose)

        # Format geo string
        self.geoStr = '{left:f}/{right:f}/{bottom:f}/{top:f}'.format(**geo.__dict__)

        # Format resolution string
        self.resStr = '{:f}'.format(geo.dx)

        # Report if requested
        if self.verbose == True:
            print('--> Geo string: {:s}'.format(self.geoStr))
            print('--> Res string: {:s}'.format(self.resStr))


    def __format_date_time__(self, date, time):
        '''
        Format date into YYYY-MM-DD and time into HH:MM:SS format.
        '''
        if self.verbose == True: print('Formatting date and time')

        # Format date
        date = '{:s}-{:s}-{:s}'.format(date[0:4], date[4:6], date[6:8])

        # Format time
        time = '{:s}:{:s}:{:s}'.format(time[0:2], time[2:4], time[4:6])

        # Format datetime string
        self.datetimeStr = '{:s}T{:s}'.format(date, time)

        # Report if requested
        if self.verbose == True:
            print('date: {:s}\ntime: {:s}'.format(date, time))
            print('--> Datetime string: {:s}'.format(self.datetimeStr))


    def __format_outname__(self, outName):
        '''
        Convert output name to GMT format.
        '''
        if self.verbose == True: print('Formatting outname')

        # Format outname
        self.outNameStr = '{:s}_%s.grd'.format(outName)

        # Report if requested
        if self.verbose == True:
            print('--> Output string: {:s}'.format(self.outNameStr))


    def __create_tide_maps__(self):
        '''
        Leverage GMT 6 to create solid Earth tide maps.
        '''
        if self.verbose == True: print('Creating solid Earth tide maps')
        
        # Format command string
        commandStr = 'gmt earthtide -T{datetimeStr:s} -G{outNameStr:s} -R{geoStr:s} -I{resStr:s} -Ce,n,v'.\
            format(**self.__dict__)

        # Report if requested
        if self.verbose == True: print('full command:\n{:s}'.format(commandStr))

        # Run command
        os.system(commandStr)

        # Report completion if requested
        if self.verbose == True: print('Tide maps generated.')


    def plot(self):
        '''
        Plot output tide maps.
        '''
        # Format tide map names
        etide = self.outNameStr.replace('%s', 'e')
        ntide = self.outNameStr.replace('%s', 'n')
        vtide = self.outNameStr.replace('%s', 'v')

        # Load data sets
        datasets = load_gdal_datasets([etide, ntide, vtide])

        # Plot data sets
        raster_multiplot(datasets.values(), ncols=3,
            cbarOrient='horizontal',
            titles=['east', 'north', 'vertical'])



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Tide maps
    tideMap = create_tide_map(dsName=inps.dsName,
        date=inps.date, time=inps.time,
        outName=inps.outName,
        verbose=inps.verbose)

    # Plot if requested
    if inps.plot == True:
        tideMap.plot()

        plt.show()