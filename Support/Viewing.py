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
from scipy.interpolate import interp1d
from osgeo import gdal
from GeoFormatting import DS_to_extent


### STATISTICS ---
def image_stats(img, mask=None, verbose=False):
    '''
    Compute the essential statistics of an image.
    '''
    # Setup
    stats = {}  # empty dictionary
    stats['Yshape'], stats['Xshape'] = img.shape

    # Mask image if mask is provided
    if mask is not None:
        img = np.ma.array(img, mask=(mask==0))
        img = img.compressed().flatten()

    # Compute statistics
    stats['mean'] = np.mean(img)
    stats['median'] = np.median(img)
    stats['min'] = np.min(img)
    stats['max'] = np.max(img)

    # Report if requested
    if verbose == True:
        print('''Image statistics
    shape: {Yshape:d} x {Xshape:d}
    mean: {mean:.3e}
    median: {median:.3e}
    min: {min:.3e}
    max: {max:.3e}'''.format(**stats))

    return stats


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


def image_clip_values(img, vmin, vmax, minPct, maxPct, verbose=False):
    '''
    Determine vmin and vmax for an image given clip values as percentiles or
     values.
    '''
    # Determine clip values
    if minPct is not None or maxPct is not None:
        minClip, maxClip = image_percentiles(img, minPct, maxPct)

    # Set clip values
    if vmin is None and minPct is not None:
        vmin = minClip  # set min clip value
    if vmax is None and maxPct is not None:
        vmax = maxClip  # set max clip value

    # Report if requested
    if verbose == True: print('Clipping image to {:f} and {:f} '.format(vmin, vmax))    

    return vmin, vmax


def image_histogram(img, mask=None, bins=128):
    '''
    Compute and show the histogram of image values.
    '''
    # Apply mask if provided
    if mask is not None:
        img = np.ma.array(img, mask=(mask==0))
        img = img.compressed().flatten()

    # Compute histogram
    hvals, hedges = np.histogram(img, bins=bins)
    hcenters = (hedges[:-1]+hedges[1:])/2

    return hcenters, hvals


def plot_histogram(img, mask=None, bins=128):
    '''
    Plot the histogram generated by image_histogram.
    '''
    # Compute histogram
    hcenters, hvals = image_histogram(img, mask, bins)

    # Spawn fresh figure
    fig, ax = plt.subplots()

    # Plot histogram
    markerline, stemlines, baseline = plt.stem(hcenters, hvals,
        linefmt='r', markerfmt='', use_line_collection=True)
    stemlines.set_linewidths(None)
    baseline.set_linewidth(0)
    ax.plot(hcenters, hvals, 'k', linewidth=2)

    # Format histogram
    ax.set_yticks([])

    return fig, ax



### SCALING ---
def equalize_image(img):
    '''
    Equalize the color balance of an image based on the histogram.
    '''
    # Compute histogram
    hcenters, hvals = image_histogram(img)
    hcenters[0], hcenters[-1] = (img.min(), img.max())

    # Integrate to build transform
    Hvals = np.cumsum(hvals)
    Hvals = Hvals - Hvals.min()  # set min value to 0
    Hvals = Hvals/Hvals.max()  # set max value to 1
    vmin, vmax = (0, 1)

    # Inverse interpolation
    Intp = interp1d(hcenters, Hvals, bounds_error=False, fill_value=0)

    img = Intp(img)

    return img, vmin, vmax



### PLOTTING ---

# >>>>>>
def plot_raster(img, mask=None, extent=None,
        cmap='viridis', cbarOrient=None,
        vmin=None, vmax=None, minPct=None, maxPct=None,
        equalize=False,
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

    # Equalize if requested
    if equalize == True:
        img, vmin, vmax = equalize_image(img)

    # Apply mask if provided
    if mask is not None:
        img = np.ma.array(img, mask=(mask==0))

    # Determine clipping values
    vmin, vmax = image_clip_values(img, vmin, vmax, minPct, maxPct)

    # Plot image
    cax = ax.imshow(img, extent=extent,
        cmap=cmap, vmin=vmin, vmax=vmax)

    # Colorbar
    if cbarOrient is not None and equalize is False:
        fig.colorbar(cax, ax=ax, orientation=cbarOrient)

    return fig, ax
# <<<<<<


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



### PROFILES ---
def plot_profile(profGeom, fig=None, ax=None):
    '''
    Plot a map profile based on the corners provided by the "profile_geometry"
     class.
    '''
    # Spawn figure if necessary
    if fig is None and ax is None:
        fig, ax = plt.subplots()

    # Plot profile
    ax.fill(profGeom.corners[:,0], profGeom.corners[:,1],
        facecolor=(0.5, 0.5, 0.5), edgecolor='k', alpha=0.5)

    # Format axis
    ax.set_aspect(1)

    return fig, ax