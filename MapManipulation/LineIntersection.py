#!/usr/bin/env python3
'''
Determine the point of intersection of two lines.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Proj


### PARSER ---
def createParser():
    Description = '''Determine the point of intersection of two lines.

This script takes input points defining the lines\' end points, then projects
 the lines into a Cartesian (UTM) coordinate system, then computes the
 intersection algebraically.
'''

    parser=argparse.ArgumentParser(description = Description,
        formatter_class = argparse.RawTextHelpFormatter)

    # Line arguments
    lineArgs = parser.add_argument_group('Line arguments')
    lineArgs.add_argument('-l1start','--line1-start', dest='l10', type=str, required=True,
        help='Line 1 start \"lon,lat\"')
    lineArgs.add_argument('-l1end', '--line1-end', dest='l11', type=str, required=True,
        help='Line 1 end \"lon,lat\"')
    lineArgs.add_argument('-l2start','--line2-start', dest='l20', type=str, required=True,
        help='Line 2 start \"lon,lat\"')
    lineArgs.add_argument('-l2end', '--line2-end', dest='l21', type=str, required=True,
        help='Line 2 end \"lon,lat\"')

    transformArgs = parser.add_argument_group('Transform arguments')
    transformArgs.add_argument('-r','--reproject', dest='reproject', default=None,
        help='''Reproject into Cartesian coordinates for better accuracy.
Default = None (don\'t reproject')
Specify Zone for UTM reprojection, or \'auto\' for automatic zone detection.''')

    outputArgs = parser.add_argument_group('Output arguments')
    outputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot outputs')

    return parser

def cmdParser(iargs = None):
    parser = createParser()

    return parser.parse_args(args = iargs)



### INTERSECTION OBJECT ---
class lineIntersect:
    '''
    Find the intersection between two lines by projecting those lines into a Cartesian coordinate system and solving
     algebraically for the point of intersection.
    '''
    def __init__(self, l10, l11, l20, l21, reproject=None):
        '''
        Solve for lines' intersection.
        '''
        # Create dictionary of line inputs
        self.linePts = {}
        self.linePts['l10'] = l10
        self.linePts['l11'] = l11
        self.linePts['l20'] = l20
        self.linePts['l21'] = l21

        # Format lines input data
        self.__checkLineInputs__()

        # Project into Cartesian coordinate system if specified
        self.projSys = reproject
        if reproject is not None:
            self.__toCartesian__(reproject)

        # Find point of intersection
        self.__findIntersection__()

    def __checkLineInputs__(self):
        '''
        Make sure line end points are in the correct (float) format.
        '''
        for key in self.linePts.keys():
            self.linePts[key] = [float(coord) for coord in self.linePts[key].split(',')]

    def __toCartesian__(self, UTMzone):
        '''
        Convert line end points to Cartesian coordinate system.
        '''
        # Determine UTM zone if 'auto'
        if UTMzone == 'auto':
            lons = [self.linePts[key][0] for key in self.linePts.keys()]
            aveLon = np.mean(lons)
            UTMzone = (int(1+(aveLon+180.0)/6.0))
        else:
            try:
                UTMzone = int(UTMzone)
            except:
                print('UTMzone not recognized')
                exit()

        print('Converting to UTM zone: {:d}'.format(UTMzone))

        # Reproject points
        self.transform = Proj(proj='utm', zone=45, ellps='WGS84')

        for key in self.linePts.keys():
            self.linePts[key] = self.transform(*self.linePts[key])

    def __findIntersection__(self):
        '''
        Find the point of intersection of the two (reprojected) lines.
        '''
        # Find equation of a line between each pair of points
        m1, b1 = self.__lineSolve__(*self.linePts['l10'], *self.linePts['l11'])
        m2, b2 = self.__lineSolve__(*self.linePts['l20'], *self.linePts['l21'])

        # Find point of intersection
        # For two lines
        #  y = m1.x + b1
        #  y = m2.x + b2
        #  x = (b2 - b1)/(m1 - m2)
        self.xIct = (b2 - b1)/(m1 - m2)
        self.yIct = m1*self.xIct+b1

        # Report intersection point in original coordinates
        if hasattr(self,'transform'):
            self.xIct, self.yIct = self.transform(self.xIct, self.yIct, inverse=True)
        print('Intersection point: {:.6f},{:.6f}'.format(self.xIct, self.yIct))

    def __lineSolve__(self, x1, y1, x2, y2):
        '''
        Find the equation of a line between points (x1, y1) and (x2, y2) of the form 
        mx + b = y, or
        beta1 x + beta2 = y.        
        '''
        m = (y2-y1)/(x2-x1)
        b = y1 - m*x1

        return m, b

    def plot(self):
        '''
        Plot input lines and results.
        '''
        # Spawn figure
        fig, ax = plt.subplots()

        # Plot lines in geo coords
        ax.plot([self.linePts['l10'][0], self.linePts['l11'][0]],
            [self.linePts['l10'][1], self.linePts['l11'][1]],
            '-ro')
        ax.plot([self.linePts['l20'][0], self.linePts['l21'][0]],
            [self.linePts['l20'][1], self.linePts['l21'][1]],
            '-bo')

        # Plot intersection point in Cartesian coordinates
        if hasattr(self,'transform'):
            self.xIct, self.yIct = self.transform(self.xIct, self.yIct)
        ax.plot(self.xIct, self.yIct, 'kd')

        # Format plot
        ax.set_aspect(1)



### MAIN ---
if __name__ == "__main__":
    # Gather input arguments
    inps = cmdParser()

    # Instantiate lineIntersect object
    L = lineIntersect(inps.l10, inps.l11, inps.l20, inps.l21,
        reproject = inps.reproject)

    # Print outputs
    if inps.plot == True:
        L.plot()

        plt.show()