#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Find and list the unique dates in a file based on various formats.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse


### PARSER ---
Description = '''Project three components of motion into satellite line of sight (LOS).
These routines work for ARIA and ISCE conventions.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')

    InputArgs.add_argument(dest='fldr', type=str,
        help='Folder in which subfolders or filenames are stored.')

    InputArgs.add_argument('-fmt','--format', dest='imgFiles', type=str, default=None,
        help='File format. Provide a file type (\'aria netcdf\', \'aria unwrappedPhase\', [None])')

    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')

    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### NAMING IDENTIFICATION ---
schemes = ['aria netcdf', 'aria unwrappedPhase']

def determine_naming_scheme(fnames, fmt=None):
    '''
    Determine the scheme for naming the files in the folder.
    '''
    # Use format if specified
    if fmt is not None:
        assert fmt in schemes, 'Specified format {:s} not is not a valid scheme ({})'.format(fmt, schemes)
        scheme = fmt
    else:
        scheme = None

    return scheme



### MAIN ---
if __name__ == '__main__':
    ## Setup
    # Gather arguments
    inps = cmdParser()

    
    ## Determine naming scheme
    determine_naming_scheme()


    print('You\'re good.')