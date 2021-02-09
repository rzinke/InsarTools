'''
SHORT DESCRIPTION
Determine a single, binary mask given an image and one or more mask values.

INHERITANCES

TESTING STATUS
Tested for background value application. Other mask types need further testing.
'''

### IMPORT MODULES ---
import os
import numpy as np
from osgeo import gdal
from scipy.stats import mode



### DETECT BACKGROUND ---
def detect_background(img, verbose=False):
    '''
    Detect background value as the most common value around the edges of an image.
    '''
    # Edge values
    edges = np.concatenate([img[:,0], img[0,:], img[:,-1], img[-1,:]])

    # Mode of edge values
    modeValues, counts = mode(edges)

    # Background value
    background = modeValues[0]

    # Report if requested
    if verbose == True: print('Background value: {:f}'.format(background))

    return background



### GENERATE MASK ---
def create_mask(img, maskArgs, verbose=False):
    '''
    Provide an image array, and one or more files or values to mask.
    Mask values come from the command line as string arguments.
    '''
    # Setup
    M, N = img.shape  # mask size
    mask = np.ones((M, N))  # everything passes preliminary mask


    if maskArgs is not None:
        # Loop through mask arguments
        for maskArg in maskArgs:
            # Determine argument type, i.e., file, value, or rule.
            if os.path.exists(maskArg):
                '''
                Assume map file.
                '''
                # Open mask image
                mskDS = open(maskArg, gdal.GA_ReadOnly)
                mskImg = mskDS.GetRasterBand(1).ReadAsArray()

                # Apply to mask
                mask *= mskImg

                # Record inferred argument type
                argType = 'map dataset'


            elif maskArg in ['bg', 'background']:
                '''
                Detect and mask background value.
                '''
                # Detect background value
                background = detect_background(img)

                # Apply to mask
                mask[img == background] = 0

                # Record inferred argument type
                argType = 'background value'


            elif maskArg[0] in ['<', '>']:
                '''
                Assume greater than/less than condition.
                '''
                # Format string to include image reference
                maskArg = 'img {:s}'.format(maskArg)

                # Evaluate
                mask[eval(maskArg)] = 0

                # Record inferred argument type
                argType = 'lt/gt condition'


            else:
                '''
                Assume single value.
                '''
                try:
                    # Try converting to float value
                    maskArg = float(maskArg)

                    # Apply to mask
                    mask[img==maskArg] = 0

                    # Record inferred argument type
                    argType = 'single value'

                except:
                    # Throw hands up \_:(_/
                    print('Mask value {} not recognized.'.format(maskArg))
                    exit()

            # Report if requested
            if verbose == True: print('Mask value {} interpreted as {:s}'.format(maskArg, argType))

    return mask