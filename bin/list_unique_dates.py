#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Find and list the unique dates in a file based on various formats.

FUTURE IMPROVEMENTS

TESTING STATUS
In development.
'''

### IMPORT MODULES ---
import argparse
import os
from glob import glob


### PARSER ---
Description = '''Project three components of motion into satellite line of sight (LOS).
These routines work for ARIA and ISCE conventions.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')

    InputArgs.add_argument(dest='srchStr', type=str,
        help='Search string leading to the subfolders or filenames are stored. E.g., \"T41/unwrappedPhase\"')

    InputArgs.add_argument('-s','--scheme', dest='scheme', type=str, default=None,
        help='File format. Provide a file type (\'aria netcdf\', \'aria unwrappedPhase\', [None])')

    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    OutputArgs.add_argument('--show-files', dest='showFiles', action='store_true',
        help='Print all detected file names.')
    OutputArgs.add_argument('--show-dates', dest='showDates', action='store_true',
        help='Show list of dates.')

    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### LOADING ---
def search_files(srchStr, verbose=False, showFiles=False):
    '''
    Determine the filenames or subfolders within the specified folder.
    '''
    # Find all results
    fpaths = (glob(srchStr))

    # Report if requested
    if verbose == True: print('{:d} files detected'.format(len(fpaths)))

    # Show all files if requested
    if showFiles == True:
        print('Files:')
        [print('\t{:s}'.format(os.path.basename(fpath))) for fpath in fpaths]

    return fpaths



### NAMING IDENTIFICATION ---
schemes = ['aria netcdf', 'aria unwrappedphase']

def determine_naming_scheme(scheme, fnames, verbose=False):
    '''
    Determine the scheme for naming the files in the folder.
    '''
    # Use format if specified
    if scheme is not None:
        scheme = scheme.lower()
        assert scheme in schemes, \
            'Specified format {:s} not is not a valid scheme ({})'.format(scheme, schemes)

    else:
        print('Format not recognized. Auto-detection not ready yet')
        exit()

    return scheme



### FORMAT DATES ---
def dates_from_fnames(scheme, fnames, verbose=False, showDates=False):
    '''
    Get a list of all dates from the file names.
    '''
    print('Formatting dates with scheme {:s}'.format(scheme))

    # Proceed according to scheme
    if scheme == 'aria unwrappedphase':
        fnames = get_basename(fnames)
        fnames = remove_extension(fnames)
        fnames = separate_underscore(fnames)
        dates = return_unique(fnames)
        dates.sort()

    else:
        print('Other schemes not supported yet.')
        exit()

    # Report if requested
    if showDates == True:
        [print(date) for date in dates]

    if verbose == True:
        print('{:d} unique dates detected'.format(len(dates)))


def get_basename(fnames):
    '''
    Return the file basename.
    '''
    fnames = [os.path.basename(fname) for fname in fnames]

    return fnames


def remove_extension(fnames):
    '''
    Remove file name extension.
    '''
    fnames = ['.'.join(fname.split('.')[:-1]) for fname in fnames]

    return fnames


def separate_underscore(fnames):
    '''
    Separate at underscore.
    '''
    elements = []
    [elements.extend(fname.split('_')) for fname in fnames]

    return elements


def return_unique(L):
    '''
    Return the unique items in a list.
    '''
    unique_items = list(set(L))

    return unique_items



### MAIN ---
if __name__ == '__main__':
    ## Setup
    # Gather arguments
    inps = cmdParser()


    ## Retrieve file names
    fnames = search_files(inps.srchStr, verbose=inps.verbose, showFiles=inps.showFiles)

    
    ## Determine naming scheme
    scheme = determine_naming_scheme(inps.scheme, fnames, verbose=inps.verbose)


    ## Gather and format dates
    # Get all dates
    dates = dates_from_fnames(scheme, fnames, verbose=inps.verbose, showDates=inps.showDates)


    print('__In development__')
