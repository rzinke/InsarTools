#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Display a profile.

FUTURE IMPROVEMENTS

TESTING STATUS
In development.
'''

### IMPORT MODULES ---
import argparse
import matplotlib.pyplot as plt
from IOsupport import load_profile_data
from Fitting import fit_atan
from Profiling import profile_binning
from Viewing import plot_profile


### PARSER ---
Description = '''Display a profile.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    InputArgs = parser.add_argument_group('INPUTS')
    InputArgs.add_argument(dest='profName', type=str,
        help='Profile name.')
    InputArgs.add_argument('--fit-type', dest='fitType', type=str, default=None,
        help='Fit function ([None], atan, linearN).')
    InputArgs.add_argument('--binning', dest='binning', action='store_true',
        help='Binning.')
    InputArgs.add_argument('--bin-width', dest='binWidth', type=float, default=None,
        help='Binning width.')
    InputArgs.add_argument('--bin-spacing', dest='binSpacing', type=float, default=None,
        help='Binning width.')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')

    return parser


def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Load profile
    profDist, profPts = load_profile_data(inps.profName, verbose=inps.verbose)


    ## Plot data
    # Spawn figure
    fig, ax = plt.subplots()

    # Determine color
    if inps.fitType is None and inps.binning is False:
        color = 'k'
    else:
        color = (0.5, 0.5, 0.5)

    # Plot profile
    ax.scatter(profDist, profPts, s=25, c=[color], label='data')

    # Binning
    if inps.binning == True:
        x, y = profile_binning(profDist, profPts,
            binWidth=inps.binWidth, binSpacing=inps.binSpacing)
        ax.plot(x, y, color='b', linewidth=3,
            label='binning')

    # Fit profile if requested
    if inps.fitType is not None:
        y, B = fit_atan(profDist, profPts, verbose=inps.verbose)
        ax.plot(profDist, y, color='r', linewidth=3,
            label='{:s} fit'.format(inps.fitType))

    # Format plot
    ax.set_xlabel('distance from start')
    ax.legend()


    plt.show()