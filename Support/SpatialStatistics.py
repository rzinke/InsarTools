'''
SHORT DESCRIPTION
Compute image spatial statistics.

FUTURE IMPROVEMENTS

TESTING STATUS
In development.
'''

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from Viewing import plot_raster


### VARIOGRAM ---
def semivariogram(img, nSamples, rSearch, nLags, mask=None, verbose=False):
    '''
    Compute the semivariogram of an image data set.
    (Semi)variance will be computed at d distances from n sample points.

    INPUTS
        img is the image to be evaluated
        mask is an image mask (0s are masked)
        nSamples is the number of samples to collect, typically a large number
        rSearch is the search radius in pixels
    '''
    if verbose == True: print('Computing semivariance')

    # Ensure nb lags < search radius in pixels
    assert nLags < rSearch, \
        'Number of lags must be less than the search radius: {:d}'.format(nSearch)

    # Parameters
    M, N = img.shape  # image dimensions

    # Query point polar coordinates
    dists = np.linspace(rSearch/nLags, rSearch, nLags)  # sample radii
    angles = np.arange(0, 2*np.pi, np.pi/15)  # angles about which to sample

    # Reshape arrays for outer product-like broadcasting
    dists = dists.reshape(-1,1)
    angles = angles.reshape(1,-1)

    # Query point arrays - preallocate coordinates for each sample point
    px = np.dot(dists, np.cos(angles)).flatten()
    py = np.dot(dists, np.sin(angles)).flatten()

    # Select n random query points in image
    Sx = np.random.randint(0, N, nSamples)
    Sy = np.random.randint(0, M, nSamples)

    if verbose == True: print('{:d} samples drawn'.format(nSamples))

    fig, ax = plt.subplots(figsize=(8,8))
    fig, ax = plot_raster(img, mask=mask, minPct=1, maxPct=99, fig=fig, ax=ax)
    ax.plot(Sx, Sy, linewidth=0, marker='.', color=(0.6,0.6,0.6))

    # First pass remove mask values
    if mask is not None:
        # Restrict to only non-masked samples
        w = (mask[Sy,Sx] != 0)  # non-masked sample indices
        Sx = Sx[w]  # valid sample x-locations
        Sy = Sy[w]  # valid sample y-locations

        # Update number of samples
        nSamples = sum(w)

        if verbose == True:
            print('{:d} samples remaining after pairing by mask'.format(nSamples))
    ax.plot(Sx, Sy, 'b.')

    # Empty arrays
    d = []
    var = []
    cov = []

    # Determine variance
    for s in range(nSamples):
        # Sample points
        sx = Sx[s]
        sy = Sy[s]

        # First pass at query points
        qx = (px + Sx[s]).astype(int)
        qy = (py + Sy[s]).astype(int)

        # Remove query points outside image bounds
        w = (qx > 0) & (qy > 0) & (qx < N) & (qy < M)
        qx = qx[w]
        qy = qy[w]

        # Remove masked query points
        w = (mask[qy,qx] != 0)
        qx = qx[w]
        qy = qy[w]

        # Recompute distances
        d = np.sqrt((qx-Sx[s])**2 + (qy-Sy[s])**2)

        # Retrieve values at query points
        qvals = img[qy,qx]

        # Calculate variances
        var = (img[sy,sx]-qvals)**2  # check this

        # Calculate covariances
        cov = img[sy,sx]*qvals

        center = (Sx[s], Sy[s])
        ax.add_patch(Circle(center, rSearch, facecolor=(0.6,0.6,0.6), alpha=0.3, zorder=2))
        ax.scatter(qx, qy, 5, color='y')