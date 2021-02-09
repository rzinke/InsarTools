'''
SHORT DESCRIPTION
Viewing functions.

INHERITANCES

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt


### MISC ---
def image_percentiles(img, minPct=1, maxPct=99, verbose=False):
    '''
    Find vmin and vmax for an image based on percentiles.
    '''
    # Determine if masked array
    if type(img) == np.ma.array:
        img = img.compressed()

    # Flatten to 1D array
    img = img.flatten()

    # Compute percentiles
    vmin, vmax = np.percentile(img, (minPct, maxPct))

    # Report if requested
    if verbose == True: print('Clipping image to {:.1f} and {:.1f} percentiles'.format(pctMin, pctMax))

    return vmin, vmax