#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Normalize a LOS velocity map to the sine of the incidence angle.

INHERITANCES
IOsupport: load_gdal_dataset, save_gdal_dataset
Masking: create_mask
GeoFormatting: DS_to_extent
Viewing: image_percentiles

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from Masking import create_mask
from GeoFormatting import DS_to_extent
from IOsupport import load_gdal_dataset, save_gdal_dataset
from Viewing import image_percentiles


### PARSER ---
Description = '''Stitch two maps belonging to two spatial data sets.'''

Examples = ''''''


def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    parser.add_argument('-v','--velocity', dest='vName', type=str,
        help='Velocity map filename.')
    parser.add_argument('-i','--incidence', dest='iName', type=str,
        help='Incidence map filename.')
    parser.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values.')
    parser.add_argument('--to-LOS', dest='normalize', action='store_false',
        help='Project back into LOS (multiply by sin(incidence angle)).')

    parser.add_argument('-o','--outName', dest='outName', type=str, default='Out', 
        help='Output name.')
    parser.add_argument('--verbose', dest='verbose', action='store_true', 
        help='Verbose mode.')
    parser.add_argument('-p','--plot', dest='plot', action='store_true', 
        help='Plot inputs and outputs.')
    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### NORMALIZATION ---
def normalize_velocity(vel, inc, mask, verbose=False):
    '''
    Normalize the velocity by the sine of the incidence angle.
    '''
    if verbose == True: print('Normalizing to incidence angle')

    # Mask to avoid dividing by zero
    vel = np.ma.array(vel, mask=(mask==0))
    inv = np.ma.array(inc, mask=(mask==0))

    # Compute sine of incidence angle
    sinInc = np.sin(np.deg2rad(inc))

    # Divide velocity by incidence angle
    normVelocity = vel/sinInc

    # Unmask
    normVelocity.mask = np.ma.nomask

    return normVelocity


def incnorm_to_LOS(vel, inc, mask, verbose=False):
    '''
    Reproject a sine-normalized map back into LOS using the incidence angles.
    '''
    if verbose == True: print('Projecting to LOS')

    # Mask rasters
    vel = np.ma.array(vel, mask=(mask==0))
    inc = np.ma.array(inc, mask=(mask==0))

    # Take sine of incidence angle
    sinInc = np.sin(np.deg2rad(inc))

    # Multiply velocity by incidnce angle
    LOSvelocity = vel*sinInc

    # Unmask
    LOSvelocity.mask = np.ma.nomask

    return LOSvelocity



### PLOTTING ---
def plot_datasets(inVelocity, incidence, outVelocity, mask, extent):
    '''
    Plot original and normalized/projected data sets.
    '''
    # Mask maps
    inVelocity = np.ma.array(inVelocity, mask=(mask==0))
    incidence = np.ma.array(incidence, mask=(mask==0))
    outVelocity = np.ma.array(outVelocity, mask=(mask==0))

    # Spawn figure
    fig, [axVel, axInc, axNorm] = plt.subplots(ncols=3)
    cbarOrient = 'horizontal'

    # Plot input velocity
    vmin, vmax = image_percentiles(inVelocity)
    caxVel = axVel.imshow(inVelocity, extent=extent,
        cmap='jet', vmin=vmin, vmax=vmax)
    axVel.set_title('Orig. veloc.')
    fig.colorbar(caxVel, ax=axVel, orientation=cbarOrient)

    # Plot incidence field
    caxInc = axInc.imshow(incidence, extent=extent)
    axInc.set_title('Incid.')
    fig.colorbar(caxInc, ax=axInc, orientation=cbarOrient)

    # Plot normalized velocity field
    caxNorm = axNorm.imshow(outVelocity, extent=extent,
        cmap='jet', vmin=vmin, vmax=vmax)
    axNorm.set_title('Normd. veloc.')
    fig.colorbar(caxNorm, ax=axNorm, orientation=cbarOrient)    



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Load data sets
    # Load velocity map
    DSvel = load_gdal_dataset(inps.vName)
    vel = DSvel.GetRasterBand(1).ReadAsArray()

    # Load incidence map
    DSinc = load_gdal_dataset(inps.iName)
    inc = DSinc.GetRasterBand(1).ReadAsArray()

    # Check velocity and incidence maps are same size
    assert vel.shape == inc.shape, 'Velocity and incidence maps must be same size'


    ## Masking
    mask = create_mask(vel, inps.maskArgs, verbose=inps.verbose)


    ## Normalize velocity map
    if inps.normalize == True:
        # Normalize
        projVel = normalize_velocity(vel, inc, mask, verbose=inps.verbose)

    elif inps.normalize == False:
        # Undo normalization (reproject back into LOS)
        projVel = incnorm_to_LOS(vel, inc, mask, verbose=inps.verbose)


    ## Save to file
    # Check filename
    if inps.outName[-4:] != '.tif': inps.outName += '.tif'

    # Save to GDAL data set
    save_gdal_dataset(inps.outName, projVel, DSvel, verbose=inps.verbose)


    ## Plot
    if inps.plot == True: plot_datasets(vel, inc, projVel, mask, DS_to_extent(DSvel))


    plt.show()