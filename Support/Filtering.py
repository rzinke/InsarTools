'''
SHORT DESCRIPTION
Filter 1D and 2D data.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import numpy as np
from scipy.signal import fftconvolve



### 2D FILTERS ---
def mean_kernel_2d(w, verbose=False):
    '''
    Apply a boxcar filter.
    INPUTS
        w is the kernel size in pixels.
    '''
    if verbose == True: print('Filtering with {:d} x {:d} mean'.format(w, w))

    # Create kernel
    h = np.ones((w,w))/w**2

    return h



### 2D FILTER APPLICATION ---
def filter_image(filtType, img, w, verbose=False):
    '''
    Apply a linear image filter via convolution.
    '''
    # Select filter
    if filtType in ['mean', 'image_mean', 'square']:
        h = mean_kernel_2d(w, verbose=verbose)

    else:
        'Filter type {:s} not supported.'.format(filtType)
        exit()

    # Apply filter
    fimg = fftconvolve(img, h, 'same')

    return fimg