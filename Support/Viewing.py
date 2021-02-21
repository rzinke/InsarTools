'''
SHORT DESCRIPTION
Viewing functions.

FUTUTRE IMPROVEMENTS
    * raster_multiplot

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
from GeoFormatting import DS_to_extent


### PLOTTING ---
def raster_multiplot(imgs, mrows=1, ncols=1,
        mask=None, extent=None,
        cmap='viridis', cbarOrient=None,
        vmin=None, vmax=None, minPct=None, maxPct=None,
        titles=None, suptitle=None,
        fig=None, ax=None):
    '''
    Plot multiple raster data sets.

    Inherits plot_raster.
    '''
    # Setup
    MN = mrows * ncols  # figure dimensions

    # Checks
    if titles is not None and len(titles) > 1:
        assert len(titles) == len(imgs), \
            'Number of titles must equal the number of images or constitute a single supertitle'

    # Spawn initial figure
    fig, axes = plt.subplots(nrows=mrows, ncols=ncols)
    if mrows > 1: axes = [ax for row in axes for ax in row]

    # Loop through images to plot
    k = 0
    for i, img in enumerate(imgs):
        # Spawn
        if (MN - k) == 0:
            # Update old figure
            fig.suptitle(suptitle)
            fig.tight_layout()

            # Spawn new figure
            fig, axes = plt.subplots(nrows=mrows, ncols=ncols)
            if mrows > 1: axes = [ax for row in axes for ax in row]
            fig.suptitle(suptitle)
            k = 0  # reset counter

        # Plot image
        fig, axes[k] = plot_raster(img, mask=mask, extent=extent,
            cmap=cmap, cbarOrient=cbarOrient,
            vmin=vmin, vmax=vmax, minPct=minPct, maxPct=maxPct,
            fig=fig, ax=axes[k])

        # Format plot
        if titles is not None: axes[k].set_title(titles[i])

        # Update counter
        k += 1

    # Format final figure
    fig.suptitle(suptitle)
    fig.tight_layout()



# >>>>>>
def plot_raster(img, mask=None, extent=None,
        cmap='viridis', cbarOrient=None,
        vmin=None, vmax=None, minPct=None, maxPct=None,
        fig=None, ax=None):
    '''
    Basic function to plot a single raster image.
    '''
    # Determine if provided image is a GDAL dataset
    if type(img) == gdal.Dataset:
        # Determine extent unless previously specified
        if extent is None:
            extent = DS_to_extent(img)

        # Retrieve image
        img = img.GetRasterBand(1).ReadAsArray()

    # Spawn figure and axis if not given
    if fig is None and ax is None:
        fig, ax = plt.subplots()

    # Mask image if mask is supplied
    if mask is not None:
        img = np.ma.array(img, mask=(mask==0))

    # Determine clip values
    if minPct is not None or maxPct is not None:
        minClip, maxClip = image_percentiles(img, minPct, maxPct)

    # Set clip values
    if vmin is None and minPct is not None:
        vmin = minClip  # set min clip value
    if vmax is None and maxPct is not None:
        vmax = maxClip  # set max clip value

    # Plot image
    cax = ax.imshow(img, extent=extent,
        cmap=cmap, vmin=vmin, vmax=vmax)

    # Colorbar
    if cbarOrient is not None:
        fig.colorbar(cax, ax=ax, orientation=cbarOrient)

    return fig, ax
# <<<<<<



### SPECIALS ---
def plot_map_mask(img, mask, fig=None, axes=None, extent=None, cmap='viridis', minPct=1, maxPct=99):
    '''
    Provide an image and a mask. Adjust parameters as desired.
    '''
    # Spawn figure and axis if necessary
    if fig is None and axes is None:
        fig, axes = plt.subplots(ncols=2)
    caxes = []  # empty color axes list
    cbarOrient='horizontal'

    # Mask image
    img = np.ma.array(img, mask=(mask==0))

    # Clip values
    vmin, vmax = image_percentiles(img, minPct=minPct, maxPct=maxPct)

    # Plot image
    caxes.append(axes[0].imshow(img, extent=extent,
        cmap=cmap, vmin=vmin, vmax=vmax))

    # Plot mask
    caxes.append(axes[1].imshow(mask, extent=extent,
        cmap='magma', vmin=0, vmax=1))

    # Format plot
    axes[0].set_title('Image')
    axes[1].set_title('Mask')
    [fig.colorbar(caxes[i], ax=axes[i], orientation=cbarOrient) for i in range(2)]

    return fig, axes



### VECTORS ---
def plot_look_vectors(Px, Py, Pz):
    '''
    Plot look vectors based on ARIA or ISCE convention.
    '''
    # Spawn figure
    fig, [axInc, axAz] = plt.subplots(ncols=2)

    # Horizontal component
    Ph = np.linalg.norm([Px, Py])

    # Plot incidence
    axInc.quiver(0, 0, Ph*np.sign(Px), Pz, color='k', units='xy', scale=1, zorder=2)

    # Plot incidence reference lines
    axInc.axhline(0, color=(0.2, 0.8, 0.5), zorder=1)
    axInc.axvline(0, color=[0.7]*3, linestyle='--', zorder=1)

    # Format incidence axis
    axInc.set_xlim([-1, 1])
    axInc.set_ylim([-0.1, 1])
    axInc.set_aspect(1)
    axInc.set_title('Incidence')

    # Plot azimuth
    axAz.quiver(0, 0, Px, Py, color='k', units='xy', scale=1, zorder=2)

    # Plot azimuth referene lines
    axAz.axhline(0, color=[0.7]*3, zorder=1)
    axAz.axvline(0, color=[0.7]*3, zorder=1)

    # Format azimuth axis
    axAz.set_xlim([-1, 1])
    axAz.set_ylim([-1, 1])
    axAz.set_aspect(1)
    axAz.set_title('Azimuth')

    fig.tight_layout()

    return fig, [axInc, axAz]



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