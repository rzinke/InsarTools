'''
SHORT DESCRIPTION
Fit various functions, including 1D signals, planes, etc.

FUTURE IMPROVEMENTS
    * PCA
    * K-means clustering
    * exponential

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from scipy import optimize


### DATE FORMATTING ---
def dates_to_datetimes(dates, datestr='%Y%m%d', verbose=False):
    '''
    Convert a list of date strings to datetime objects.
    '''
    if verbose == True:
        print('Converting {:d} dates to datetimes with format {:s}'.format(len(dates), datestr))

    # Convert to datetimes
    datetimes = [dt.datetime.strptime(date, datestr) for date in dates]

    return datetimes


def time_since_reference(datetimes, refDate=None, verbose=False):
    '''
    Find the time difference in days between each date and the reference time.
    '''
    # Determine reference date
    if refDate is None: refDate = datetimes[0]

    # Calculate time differences
    timedeltas = [(datetime-refDate).days/365.25 for datetime in datetimes]

    return timedeltas



### TIMESERIES FITTING ---
def fit_linear(x, y, verbose=False, plot=False):
    '''
    Fit an offset and secular trend to a 1D signal.
    '''
    # Parameters
    assert len(x) == len(y), 'Arrays x and y must be the same length'
    N = len(x)

    # Design matrix
    G = np.ones((N, 2))
    G[:,1] = x.flatten()

    # Invert for fit parameters
    B = np.linalg.inv(np.dot(G.T, G)).dot(G.T).dot(y)

    # Reconstruct curve
    yhat = G.dot(B)

    # Compute residuals
    res = y - yhat
    RSS = np.sqrt(np.sum(res**2))  # root sum of squares
    expRes = RSS/len(res)  # mean residual

    # Report if requested
    if verbose == True:
        B = B.flatten()
        print('Linear fit')
        print('x^0 {:f}\nx^1 {:f}'.format(*B))
        print('Expected residual: {:f}'.format(expRes))

    # Plot if requested
    if plot == True:
        fig, ax = plt.subplots()
        ax.plot(x, y, 'k.', label='data')
        ax.plot(x, yhat, 'b', label='fit')

    return yhat, B


def fit_periodic(x, y, freq=1, verbose=False, plot=False):
    '''
    Fit a periodic, linear, and offset to a 1D signal.
    Assumes a period/freq of 1 year.
    '''
    # Parameters
    assert len(x) == len(y), 'Arrays x and y must be the same length'
    N = len(x)

    # Design matrix
    G = np.ones((N, 4))
    G[:,1] = x.flatten()
    G[:,2] = np.sin(2*np.pi*freq*x).flatten()
    G[:,3] = np.cos(2*np.pi*freq*x).flatten()

    # Invert for fit parameters
    B = np.linalg.inv(np.dot(G.T, G)).dot(G.T).dot(y)

    # Reconstruct curve
    yhat = G.dot(B)

    # Compute residuals
    res = y - yhat
    RSS = np.sqrt(np.sum(res**2))  # root sum of squares
    expRes = RSS/len(res)  # mean residual

    # Report if requested
    if verbose == True:
        B = B.flatten()
        print('Periodic fit')
        print('x^0 {:f}\nx^1 {:f}\nsin(x) {:f}\ncos(x) {:f}'.format(*B))
        print('Expected residual: {:f}'.format(expRes))

    # Plot if requested
    if plot == True:
        fig, ax = plt.subplots()
        ax.plot(x, y, 'k.', label='data')
        ax.plot(x, yhat, 'b', label='fit')

    return yhat, B


def fit_atan(x, y, initial=None, verbose=False, plot=False):
    '''
    Fit an arctangent function to a series of data points.
    '''
    # Model function
    atan = lambda x, a, b, c, d: a*np.arctan(b*(x-c))+d

    # Fit to data points
    B, cov = optimize.curve_fit(atan, x, y, p0=initial)

    # Reconstruct curve
    yhat = atan(x, *B)

    # Compute residuals
    res = y - yhat
    RMSE = np.sqrt(np.sum(res**2)/len(x))  # root mean squared error

    # Report if requested
    if verbose == True:
        a, b, c, d = B
        a = a*np.pi/2
        print('Arctangent fit')
        print('a: {:.4f}\nb: {:.4f}\nc: {:.4f}\nd: {:.4f}'.format(a, b, c, d))
        print('RMSE: {:f}'.format(RMSE))

    # Plot if requested
    if plot == True:
        fig, ax = plt.subplots()
        ax.plot(x, y, 'k.', label='data')
        ax.plot(x, yhat, 'r', label='fit')

    return yhat, B



### PRINCIPAL COMPONENT ANALYSIS ---



### K-MEANS CLUSTER ANALYSIS ---



### SURFACE FITTING ---
def design_matrix2d(x, y, degree, verbose=False):
    '''
    Build a design matrix for a set of points that vary in x and y.
    '''
    if verbose == True: print('Building design matrix')

    # Parameters
    nx = len(x)
    ny = len(y)
    assert nx == ny, \
        'Length of x ({:d}) and y({:d}) must be the same'.format(nx, ny)
    n = nx  # number of data points

    # Design matrix - depends on polynomial degree
    G = np.ones((n, 1+2*degree))  # empty design matrix
    for i in range(1, degree+1):
        # Populate design matrix
        G[:,2*i-1] = x**i
        G[:,2*i] = y**i

    return G


def fit_surface_to_points(x, y, z, degree=1, verbose=False):
    '''
    Fit a n-degree polynomial to the given points.
    x, y, and z are 1D arrays.
    '''
    if verbose == True: print('Fitting surface to points')

    # Design matrix
    G = design_matrix2d(x, y, degree=degree, verbose=verbose)

    # Plane fit parameters
    B = np.linalg.inv(np.dot(G.T, G)).dot(G.T).dot(z)

    # Reconstruct z
    zhat = G.dot(B)

    # Compute residuals
    res = z - zhat  # residuals
    rms = np.sqrt(np.mean(res**2))

    # Report if requested
    if verbose == True:
        B = B.flatten()
        print('Polynomial fit')
        print('Order: {:d}'.format(degree))
        print('\toffset: {:f}'.format(B[0]))
        [print('\tX^{i:d} {X:f} Y^{i:d} {Y:f}'.format(**{'i':i, 'X':B[2*i-1], 'Y':B[2*i]})) \
            for i in range(1, degree+1)]
        print('RMS residual: {:f}'.format(rms))

    return zhat, B


def fit_surface_to_points_weighted(x, y, z, weights, degree=1, verbose=False):
    '''
    Fit a n-degree polynomial to the given points.
    x, y, z, and weights are 1D arrays.
    '''
    if verbose == True: print('Fitting surface to points (with weighting)')

    # Check inputs
    assert len(x) == len(y) == len(z) == len(weights), 'Input lengths not the same'
    n = len(z)  # number of data points

    # Design matrix
    G = design_matrix2d(x, y, degree=degree, verbose=verbose)

    # Weighting matrix
    W = np.identity(n)
    for i in range(n):
        W[i,i] = weights[i]

    # Plane fit parameters
    B = np.linalg.inv(G.T.dot(W).dot(G)).dot(G.T).dot(W).dot(z)

    # Reconstruct z
    zhat = G.dot(B)

    # Compute residuals
    res = z - zhat  # residuals
    rms = np.sqrt(np.mean(res**2))

    # Report if requested
    if verbose == True:
        B = B.flatten()
        print('Polynomial fit')
        print('Order: {:d}'.format(degree))
        print('\toffset: {:f}'.format(B[0]))
        [print('\tX^{i:d} {X:f} Y^{i:d} {Y:f}'.format(**{'i':i, 'X':B[2*i-1], 'Y':B[2*i]})) \
            for i in range(1, degree+1)]
        print('RMS residual: {:f}'.format(rms))

    return zhat, B


def fit_surface_to_image(img, mask, degree=1, dx=1, dy=1, decimation=0, verbose=False):
    '''
    Fit a n-degree polynomial to the image data set.
    Exclude mask values. Decimate by 10^[decimation].
    Horizontal units are arbitrary as long as pixels dimensions are equal.
    '''
    # Decimation factor
    dFactor = int(10**decimation)

    # Image parameters
    M, N = img.shape
    MN = M*N

    # Reformat image as MN x 1 array
    Z = img.copy()
    Z = Z.reshape(MN, 1)

    # Establish grid
    y = np.linspace(-1, 1, M)/dy
    x = np.linspace(-1, 1, N)/dx
    X, Y = np.meshgrid(x, y)

    # Design matrix - depends on polynomial degree
    G = design_matrix2d(X.flatten(), Y.flatten(), degree=degree, verbose=verbose)

    # Reduce data to valid data points
    nds = mask.flatten()==1  # non-masked indices
    Gvld = G[nds,:]  # valid data points
    Zvld = Z[nds]  # valid data points

    # Decimate for speed
    Gvld = Gvld[::dFactor]
    Zvld = Zvld[::dFactor]

    # Plane fit parameters
    B = np.linalg.inv(np.dot(Gvld.T, Gvld)).dot(Gvld.T).dot(Zvld)

    # Reconstruct surface
    surface = G.dot(B).reshape(M, N)

    # Compute residuals
    res = img - surface  # residuals
    res = np.ma.array(res, mask=(mask==0))  # mask residuals
    res = res.compressed().flatten()  # use only valid residuals
    rms = np.sqrt(np.mean(res**2))  # root sum of squares

    # Report if requested
    if verbose == True:
        B = B.flatten()
        print('Polynomial fit')
        print('Order: {:d}'.format(degree))
        print('\toffset: {:f}'.format(B[0]))
        [print('\tX^{i:d} {X:f} Y^{i:d} {Y:f}'.format(**{'i':i, 'X':B[2*i-1], 'Y':B[2*i]})) \
            for i in range(1, degree+1)]
        print('RMS residual: {:f}'.format(rms))

    return surface, B



### VELOCITY MAP ---
class velocity_from_timeseries:
    def __init__(self, disps, times, mask=None, fitType='linear', verbose=False):
        '''
        Compute velocity maps from a displacement timeseries.

        INPUTS
            disps is a nEpochs x M x N numpy array of displacement values
            times is a nEpochs-length list of datetime time deltas
        '''
        # Parameters
        self.verbose = verbose
        self.nEpochs, self.M, self.N = disps.shape

        # Check inputs
        assert self.nEpochs == len(times), 'Number of displacement epochs must equal number of times.'

        # Check fit type and determine function
        self.__format_fit_type__(fitType)

        # Build design matrix
        self.__build_design_matrix__(times)

        # Invert design matrix
        self.__invert_design_matrix__()

        # Compute velocities
        self.__compute_velocities__(disps)


    def __format_fit_type__(self, fitType):
        '''
        Check and format the type of fit to apply.
        '''
        # Lower case
        fitType = fitType.lower()

        # Check that type is recognized
        assert fitType in ['linear', 'seasonal'], \
            'Fit type {:s} is not valid, use \'linear\' or \'seasonal\''.format(fitType)

        # Record value
        self.fitType = fitType

        # List of map components
        if self.fitType == 'linear':
            self.componentsList = ['offsetMap', 'velocityMap', 'residualMap']
        elif self.fitType == 'seasonal':
            self.componentsList = ['offsetMap', 'velocityMap', 'sineAmpMap', 'sinePhaseMap', 'residualMap']


    def __build_design_matrix__(self, times):
        '''
        Create the design matrix based on the fit type.
        First component is always offset; second is always linear velocity.
        '''
        # For linear fit
        if self.fitType == 'linear':
            # Formulate design matrix
            self.G = np.ones((self.nEpochs, 2))
            self.G[:,1] = times.flatten()

        # For seasonal fit
        elif self.fitType == 'seasonal':
            # Frequency
            freq = 1  # per year

            # Formulate design matrix
            self.G = np.ones((self.nEpochs, 4))
            self.G[:,1] = times.flatten()
            self.G[:,2] = np.cos(2*np.pi*freq*times).flatten()
            self.G[:,3] = np.sin(2*np.pi*freq*times).flatten()


    def __invert_design_matrix__(self):
        '''
        Invert the design matrix based on the fit type.
        Ideally only do this once.
        '''
        # Invert design matrix
        self.inverseMatrix = np.linalg.inv(np.dot(self.G.T, self.G)).dot(self.G.T)


    def __find_velocity_components__(self, displacements):
        '''
        Solve for the beta vector by taking the dot product Ginv.displacements.
        The second term in the beta vector is the linear component.

        Feed this through a loop to parallelize.
        '''
        B = self.inverseMatrix.dot(displacements)

        return B


    def __compute_velocities__(self, disps):
        '''
        Solve for velocities using the given fit type.
        '''
        if self.verbose == True:
            print('Computing velocity with fit type: {:s}'.format(self.fitType))

        # Setup
        solutions = []
        self.offsetMap    = np.zeros((self.M, self.N))
        self.velocityMap  = np.zeros((self.M, self.N))
        self.sineAmpMap   = np.zeros((self.M, self.N))
        self.sinePhaseMap = np.zeros((self.M, self.N))
        self.residualMap  = np.zeros((self.M, self.N))

        # Loop through each pixel
        for i in range(self.M):
            for j in range(self.N):
                # Recall displacements
                pxDisps = disps[:,i,j]

                # Solve for component scaling factors
                B = self.__find_velocity_components__(pxDisps)

                # Store linear solutions
                self.offsetMap[i,j] = B[0]
                self.velocityMap[i,j] = B[1]
                if self.fitType == 'seasonal':
                    # Sinusoid amplitude
                    self.sineAmpMap[i,j] = np.sqrt(B[2]**2 + B[3]**2)

                    # Sinusoid phase
                    self.sinePhaseMap[i,j] = np.arctan2(B[3], B[2])

                # Compute and store residuals
                residuals = pxDisps - self.G.dot(B)
                self.residualMap[i,j] = np.sqrt(np.mean(residuals**2))