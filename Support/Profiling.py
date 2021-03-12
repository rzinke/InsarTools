'''
SHORT DESCRIPTION
Image profile support.

FUTURE IMPROVEMENTS

TESTING STATUS
In development.

POLYLINE NAMING CONVENTION
                    profile
                       o
                       | }-spacing
                       o
                       |
                       o
                       |
                 __.___.___.__. polyline
                /      |
               .       o
              /        |
    width    /         o
    |   |   .          |
 ___.___.__/           o
 \ /       
 offset

PROFILE NAMING CONVENTION
   _____
  |     |
  |     | \
  | -w- | length
  |     | /
  |     |
--x--o--x- anchor
  |w2 w2|
  |     |
'''

### IMPORT MODULES ---
import numpy as np
from scipy.interpolate import interp1d


### VECTOR MATH ---
def rotation_matrix(theta):
    '''
    Definition of a standard rotation matrix.
    CCW is positive.
    Theta given in radians.
    '''
    # Rotation matrix
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])

    return R


def determine_pointing_vector(pxStart, pyStart, pxEnd, pyEnd, verbose=False):
    '''
    Find the pointing vector between the two points and the original vector
     length.
    '''
    # Find pointing vector
    p = np.array([pxEnd-pxStart, pyEnd-pyStart])
    pLen = np.linalg.norm(p)
    p = p/pLen

    # Report if requested
    if verbose == True:
        print('unit vector: {:.3f}, {:.3f}'.format(*p))
        print('vector length: {:f}'.format(pLen))

    return p, pLen


def pointing_vector_to_angle(px, py, verbose=False):
    '''
    Convert pointing vector px, py to angle.
    Theta returned in radians.
    '''
    theta = np.arctan2(py, px)

    return theta


def rotate_coordinates(X, Y, theta, verbose=False):
    '''
    Rotate X, Y coordinates by the angle theta.
    '''
    # Setup
    M, N = X.shape
    MN = M*N

    R = rotation_matrix(theta)

    # Reshape coordinate points into 2 x MN array for rotation
    C = np.vstack([X.reshape(1, MN),
                   Y.reshape(1, MN)])
    del X, Y  # clear memory

    # Rotate coordinates
    C = R.dot(C)

    # Reshape coordinate matrices
    X = C[0,:].reshape(M,N)
    Y = C[1,:].reshape(M,N)
    del C  # clear memory

    return X, Y



### DISTANCE ALONG PROFILE ---
def polyline_length(x, y):
    '''
    Calculate the distance along a polyline.
    '''
    # Parameters
    assert len(x) == len(y)
    N = len(x)  # number of data points
    distances = np.zeros(N)  # empty array of distances

    # Loop through each point
    for i in range(N-1):
        # Distance components
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]

        # Euclidean distance
        d = np.sqrt(dx**2 + dy**2)

        distances[i+1] = distances[i] + d  # start with zero, add cumulative

    return distances


def points_along_polyline(lineX, lineY, width, offset=0, verbose=False):
    '''
    Compute the x,y coordinates of points along a polyline at the specified
     spacing (width).

    INPUTS
        lineX, lineY are the coordinates of the polyline vertices
        width is the profile width
        offset is the distance between the first vertex and the first point

    OUTPUTS
        qx, qy are the profile points
    '''
    # Compute the distance array along the profile
    l = polyline_length(lineX, lineY)

    # Array of query points (function of distance)
    q = np.arange(offset, l[-1], width)

    # Interpolate along distance axis
    I = interp1d(l, np.column_stack([lineX, lineY]), axis=0)  # interp function

    # Solve for query point coordinates
    coords = I(q)
    qx = coords[:,0]  # x-coordinates
    qy = coords[:,1]  # y-coordinates

    return qx, qy


def find_profile_anchors(lineX, lineY, width, verbose=False):
    '''
    Find the anchor points of profiles along a polyline.
    '''
    # Find start anchors (0 offset)
    startAnchors = points_along_polyline(lineX, lineY, width, verbose=verbose)
    startAnchors = np.column_stack(startAnchors)

    # Find end anchors (w offset)
    endAnchors = points_along_polyline(lineX, lineY, width, offset=width, verbose=verbose)
    endAnchors = np.column_stack(endAnchors)

    # Clip start anchors to same number as end anchors
    startAnchors = startAnchors[:endAnchors.shape[0],:]

    return startAnchors, endAnchors


def find_profile_geometries(lineX, lineY, profWidth, profLen, verbose=False):
    '''
    Determine the coordinates of profile vertices along a polyline.
    '''
    # Setup
    profCoords = {}

    # Find the profile anchors along the polyline
    startAnchors, endAnchors = find_profile_anchors(lineX, lineY, profWidth, verbose=verbose)

    # Create profGeom objects to describe the profile coordinates
    nAnchors = len(startAnchors)
    profGeoms = [profGeom(startAnchors[i], endAnchors[i], profLen) for i in range(nAnchors)]

    return profGeoms


class profGeom:
    def __init__(self, startAnchor, endAnchor, profLen):
        '''
        Store the coordinates and properties composing a profile.
        Provide the starting and ending anchor points.
        Automatically calculate the orientation vectors.
        '''
        self.startAnchor = startAnchor
        self.endAnchor = endAnchor

        # Determine mid points
        self.__find_midpoint__()

        # Determine vectors
        self.__determine_vectors__()

        # Determine profile start and end
        self.__determine_start_end__(profLen)

        # Determine profile corners
        self.__determine_corners__(profLen)


    def __find_midpoint__(self):
        '''
        Find the midpoint of the profile.
        '''
        self.midAnchor = (self.endAnchor - self.startAnchor)/2 + self.startAnchor


    def __determine_vectors__(self):
        '''
        Determine orientation vector and orthogonal orientation.
        '''
        self.orient, _ = determine_pointing_vector(*self.startAnchor, *self.endAnchor)
        self.orthog = np.array([-self.orient[1], self.orient[0]])


    def __determine_start_end__(self, profLen):
        '''
        Determine the start and end points of a profile.
        '''
        self.profStart = profLen*self.orthog+self.midAnchor
        self.profEnd = -profLen*self.orthog+self.midAnchor


    def __determine_corners__(self, profLen):
        '''
        Determine profile corners.
        '''
        corner1 =  profLen*self.orthog+self.startAnchor
        corner2 =  profLen*self.orthog+self.endAnchor
        corner3 = -profLen*self.orthog+self.endAnchor
        corner4 = -profLen*self.orthog+self.startAnchor

        self.corners = np.vstack([corner1, corner2, corner3, corner4, corner1])



### PROFILING ---
def extract_profile(img, pxStart, pyStart, pxEnd, pyEnd, width, mask=None, verbose=False):
    '''
    Extract polyline-perpendicular profiles from one or more images.
    Width given in pixels
    '''
    if verbose == True: print('Extracting profile')

    # Parameters
    M, N = img.shape
    MN = M*N
    w2 = int(width/2)
    copyimg = img.copy()

    # Build grid
    x = np.arange(N)
    y = np.arange(M)

    X, Y = np.meshgrid(x, y)

    # Recenter at starting coordinates
    X = X - pxStart
    Y = Y - pyStart

    # Determine direction unit vector
    p, pLen = determine_pointing_vector(pxStart, pyStart, pxEnd, pyEnd, verbose=verbose)

    # Determine rotation angle
    theta = pointing_vector_to_angle(p[0], p[1], verbose=verbose)

    # Rotate coordinates
    X, Y = rotate_coordinates(X, Y, -theta)

    # Extract valid pixels
    validPts = np.ones(img.shape)
    validPts[X>pLen] = 0
    validPts[X<0] = 0
    validPts[Y<-w2] = 0
    validPts[Y>w2] = 0

    if mask is not None:
      validPts[mask == 0] = 0

    # Extract valid points
    profDist = X[validPts == 1].flatten()  # point distances
    profPts = img[validPts == 1].flatten()  # valid points
    del validPts  # clear memory

    return profDist, profPts