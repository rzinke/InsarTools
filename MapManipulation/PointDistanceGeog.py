#!/usr/bin/env python3
'''
Find the distance from one point to another
'''

### IMPORT MODULES ---
import argparse
import numpy as np
from pyproj import Geod


### PARSER ---
def createParser():
    Description = '''Given two points in GEOGRAPHIC coordinates, calculate the distance from one to another.
'''

    parser=argparse.ArgumentParser(description = Description,
        formatter_class = argparse.RawTextHelpFormatter)

    # Point arguments
    pointArgs = parser.add_argument_group('Point arguments')
    pointArgs.add_argument('-p1','--point1', dest='p1', type=str, required=True,
        help='First point \"lon,lat\"')
    pointArgs.add_argument('-p2','--point2', dest='p2', type=str, required=True,
        help='Second point \"lon,lat\"')

    return parser

def cmdParser(iargs = None):
    parser = createParser()

    return parser.parse_args(args = iargs)



### FIND DISTANCE ---
def findDistance(p1, p2):
    '''
    Find the distance between two points in lat/lon.
    '''
    # Format point inputs
    lon1, lat1 = [float(p) for p in p1.split(',')]
    lon2, lat2 = [float(p) for p in p2.split(',')]

    # Compute distance using Geod object
    g = Geod(ellps='WGS84')
    az12, az21, dist_m = g.inv(lon1, lat1, lon2, lat2)

    # Convert distance to degrees
    dist_deg = dist_m/111000

    # Report distance
    print('Distance (m): {:f}'.format(dist_m))
    print('Distance (deg): {:f}'.format(dist_deg))



### MAIN ---
if __name__ == '__main__':
    # Gather input arguments
    inps = cmdParser()

    # Find distance between two points
    findDistance(inps.p1, inps.p2)