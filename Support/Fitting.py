'''
SHORT DESCRIPTION
Fit various functions, including 1D signals, planes, etc.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from scipy import optimize


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
    a, b, c, d = B
    a = a*np.pi/2

    if verbose == True:
        print('Arctangent fit')
        print('a: {:.4f}\nb: {:.4f}\nc: {:.4f}\nd: {:.4f}'.format(a, b, c, d))
        print('RMSE: {:f}'.format(RMSE))

    # Plot if requested
    if plot == True:
        fig, ax = plt.subplots()
        ax.plot(x, y, 'k.', label='data')
        ax.plot(x, yhat, 'r', label='fit')

    return yhat, B



### SURFACE FITTING ---
def fit_surface(img, mask, degree=1, dx=1, dy=1, decimation=0, verbose=False, plot=False):
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
    G = np.ones((MN, 1+2*degree))  # empty design matrix
    for i in range(1, degree+1):
        # Populate design matrix
        G[:,2*i-1] = X.flatten()**i
        G[:,2*i] = Y.flatten()**i

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
    RSS = np.sqrt(np.sum(res**2))  # root sum of squares
    expRes = RSS/len(res)  # mean residual

    # Report if requested
    if verbose == True:
        B = B.flatten()
        print('Polynomial fit')
        print('Order: {:d}'.format(degree))
        print('\toffset: {:f}'.format(B[0]))
        [print('\tX^{i:d} {X:f} Y^{i:d} {Y:f}'.format(**{'i':i, 'X':B[2*i-1], 'Y':B[2*i]})) \
            for i in range(1, degree+1)]
        print('Expected residual: {:f}'.format(expRes))

    # Plot if requested
    if plot == True:
        fig, ax = plt.subplots()
        ax.imshow(surface)
        ax.set_title('Surface')

    return surface, B



### MISCELLANEOUS ---
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