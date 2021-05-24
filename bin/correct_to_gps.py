#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Fit a displacement or velocity field to a set of GPS measurements.

FUTURE IMPROVEMENTS

TESTING STATUS
Tested.
'''

### IMPORT MODULES ---
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import confirm_outdir, confirm_outname_ext, GPSdata, load_gdal_dataset, save_gdal_dataset
from Masking import create_mask
from GeoFormatting import get_raster_size, lola_to_xy, xy_to_lola, transform_to_extent
from RasterResampling import sample_points_from_raster
from Fitting import fit_surface_to_points_weighted, design_matrix2d
from Viewing import image_percentiles, plot_raster


### PARSER ---
Description = '''Fit a displacement or velocity field to a set of GPS measurements.'''

Examples = ''''''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    ImgArgs = parser.add_argument_group('IMAGE INPUTS')
    ImgArgs.add_argument('-i','--image-file', dest='imgName', type=str, required=True,
        help='Image file name.')
    ImgArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')
    ImgArgs.add_argument('--coherence', dest='cohName', type=str, default=None,
        help='Coherence matrix for weighting values. ([None]).')


    GPSargs = parser.add_argument_group('GPS INPUTS')
    GPSargs.add_argument('-g','--gps-file', dest='GPSname', type=str, required=True,
        help='GPS file name.')
    GPSargs.add_argument('-lon','--lon-column', dest='lonCol', type=int, default=0,
        help='Longitude column')
    GPSargs.add_argument('-lat','--lat-column', dest='latCol', type=int, default=1,
        help='Latitude column')
    GPSargs.add_argument('-c','--data-column', dest='DATAcol', type=int, default=2,
        help='Data column (starts at 0, [2]).')
    GPSargs.add_argument('--header-rows', dest='headerRows', type=int, default=0,
        help='Number of header rows to skip ([0], 1, 2, ...)')
    GPSargs.add_argument('--delimiter', dest='delimiter', type=str, default=' ',
        help='Text file delimiter ([\' \'],... etc.)')
    GPSargs.add_argument('--scale-factor', dest='scaleFactor', type=float, default=1,
        help='Scale factor for GPS, e.g., 0.001 converts mm to m. ([1]).')


    FitArgs = parser.add_argument_group('FIT INPUTS')
    FitArgs.add_argument('-d','--fit-degree', dest='fitDegree', type=int, default=3,
        help='Surface fit degree (0, [1], 2, 3, ...).')
    FitArgs.add_argument('-r','--search-radius', dest='searchRadius', type=int, default=0,
        help='Search radius around GPS station (pixels). Default 0 uses the closest single pixel.')
    FitArgs.add_argument('--gps-clip', dest='GPSclip', nargs=4, type=float, default=None,
        help='Pre-clip the GPS stations by location. (west east south north). ([None]).')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true',
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot results.')
    OutputArgs.add_argument('-o','--outname', dest='outName', type=str, default='Adjusted',
        help='GPS-corrected field.')

    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### GPS CORRECTION ---
class GPScorrection:
    def __init__(self, verbose=False):
        '''
        Correct an image to GPS point measurements.
        '''
        # Parameters
        self.verbose = verbose


    ## Loading
    def load_gps(self, GPSname, lonCol=0, latCol=1, DATAcol=2, headerRows=0, delimiter=' ', scaleFactor=1):
        '''
        Load GPS data from file.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Loading GPS data')

        # Load data using GPSdata object
        self.GPS = GPSdata(verbose=self.verbose)
        self.GPS.load_from_file(GPSname, lonCol=lonCol, latCol=latCol, DATAcol=DATAcol,
            headerRows=headerRows, delimiter=delimiter)

        # Scale by factor
        self.GPS.data = self.GPS.data*scaleFactor


    def load_image(self, imgName, cohName=None, maskArgs=None):
        '''
        Load image file.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Loading image data')

        self.__load_image_data__(imgName)

        # Load or create coherence matrix
        self.__load_coherence_matrix__(cohName)

        # Create mask
        self.mask = create_mask(self.img, maskArgs, verbose=self.verbose)


    def __load_image_data__(self, imgName):
        '''
        Load the image name and parse the spatial information.
        '''
        # Load image data set
        DS = load_gdal_dataset(imgName, verbose=self.verbose)
        self.proj = DS.GetProjection()
        self.tnsf = DS.GetGeoTransform()
        self.M, self.N = get_raster_size(DS)
        self.extent = transform_to_extent(self.tnsf, self.M, self.N, verbose=self.verbose)
        self.img = DS.GetRasterBand(1).ReadAsArray()

        del DS  # clear memory


    def __load_coherence_matrix__(self, cohName):
        '''
        Load coherence matrix or create uniform matrix if none specified.
        '''
        if cohName is not None:
            if self.verbose == True:
                print('Loading coherence data from file')

            # Load map data
            DScoh = load_gdal_dataset(cohName, verbose=self.verbose)

            # Check that map dimensions are the same
            M, N = get_raster_size(DScoh)
            assert (M == self.M) & (N == self.N), \
                'Coherence map size ({:d} x {:d}) is not the same as image size ({:d} x {:d})'.format(M, N, self.M, self.N)

            # Retrieve coherence image data
            self.coh = DScoh.GetRasterBand(1).ReadAsArray()

            # Clear memory
            del DScoh

        else:
            if self.verbose == True:
                print('Using coherence 1.0')
            self.coh = np.ones((self.M, self.N))


    ## GPS correction
    def correct_to_gps(self, GPSclip=None, fitDegree=1, searchRadius=0):
        '''
        Correct the image to the GPS data.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Correcting to GPS')

        # Clip to specified extent
        if GPSclip is not None: self.GPS.crop_to_extent(*GPSclip)

        # Clip to map bounds
        self.GPS.crop_to_extent(*self.extent)

        # Compute GPS - image difference
        self.__compute_difference__(searchRadius)

        # Create difference ramp
        self.__create_difference_ramp__(fitDegree)

        # Adjust image to GPS
        self.__adjust_image__(searchRadius)


    def __compute_difference__(self, searchRadius):
        '''
        Compute the difference between the GPS and the image data.
        '''
        if self.verbose == True: print('Computing difference (GPS - img)')

        # Convert GPS coordinates to image pixel coordinates
        px, py = lola_to_xy(self.tnsf, self.GPS.lon, self.GPS.lat, verbose=self.verbose)

        # Extract sample points from image
        self.sX, self.sY, sImg, sNdx = sample_points_from_raster(self.img, px, py, mask=self.mask,
            searchR=searchRadius, verbose=self.verbose)

        # Remove masked values from GPS data set
        self.GPS.crop_by_index(sNdx)

        # Compute differences GPS - img
        self.diffs = self.GPS.data - sImg

        # Extract weighting values
        _, _, self.weights, _ = sample_points_from_raster(self.coh, px, py, mask=self.mask,
            searchR=searchRadius, verbose=self.verbose)


    def __create_difference_ramp__(self, fitDegree):
        '''
        Creat a surface fitting the differences.
        '''
        # Fit ramp to differences
        _, rampParams = fit_surface_to_points_weighted(self.sX, self.sY, self.diffs,
            weights=self.weights,
            degree=fitDegree,
            verbose=self.verbose)

        # Create image grid
        imgx = np.arange(self.N)
        imgy = np.arange(self.M)
        imgX, imgY = np.meshgrid(imgx, imgy)

        # Create full surface based on fit parameters
        G = design_matrix2d(imgX.flatten(), imgY.flatten(), degree=fitDegree, verbose=self.verbose)
        diffRamp = G.dot(rampParams)
        self.diffRamp = diffRamp.reshape(self.M, self.N)


    def __adjust_image__(self, searchRadius):
        '''
        Adjust image by adding difference ramp to the image.
            diff = GPS - img
         => img + diff = GPS
        '''
        if self.verbose == True: print('Adjusting image')

        # Adjust image
        self.imgAdjusted = self.img + self.diffRamp

        # Compute residuals
        _, _, sAdjusted, _ = sample_points_from_raster(self.imgAdjusted, self.sX, self.sY,
            mask=self.mask, searchR=searchRadius, verbose=self.verbose)
        self.resids = self.GPS.data - sAdjusted
        RMSresid = np.sqrt(np.sum(self.resids**2))

        # Report if requested
        if self.verbose == True: print('RMS residual after fit: {:f}'.format(RMSresid))


    ## Saving
    def save(self, outName):
        '''
        Save image to file.
        '''
        if inps.verbose == True:
            print('*'*32)
            print('Saving outputs')

        # Format output name
        outName = confirm_outname_ext(outName, ext=['tif'], verbose=self.verbose)
        confirm_outdir(outName)

        # Save adjusted image to file
        save_gdal_dataset(outName, self.imgAdjusted, mask=self.mask, proj=self.proj, tnsf=self.tnsf,
            verbose=self.verbose)


    ## Plotting
    def plot(self):
        '''
        Plot inputs and outputs.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Plotting results')

        # Spawn figure and axes
        fig, [axOG, axDiff, axCorr] = plt.subplots(figsize=(8,7), ncols=3)
        cbarOrient = 'horizontal'

        # Plot inputs
        vmin, vmax = image_percentiles(self.img, minPct=1, maxPct=99, verbose=False)

        plot_raster(self.img, mask=self.mask, extent=self.extent,
            cmap='viridis', cbarOrient=cbarOrient,
            vmin=vmin, vmax=vmax,
            fig=fig, ax=axOG)

        axOG.scatter(self.GPS.lon, self.GPS.lat, 16, c=self.GPS.data,
            cmap='viridis', vmin=vmin, vmax=vmax)

        axOG.set_title('Orig data')

        # Plot differences
        vmin, vmax = image_percentiles(self.diffRamp, minPct=1, maxPct=99, verbose=False)

        plot_raster(self.diffRamp, mask=self.mask, extent=self.extent,
            cmap='jet', cbarOrient=cbarOrient,
            vmin=vmin, vmax=vmax,
            fig=fig, ax=axDiff)

        axDiff.scatter(self.GPS.lon, self.GPS.lat, 16, c=self.diffs,
            cmap='jet', vmin=vmin, vmax=vmax)

        axDiff.set_title('Differences')


        # Plot adjusted image and residuals
        plot_raster(self.imgAdjusted, mask=self.mask, extent=self.extent,
            cmap='viridis', cbarOrient=cbarOrient,
            minPct=1, maxPct=99,
            fig=fig, ax=axCorr)

        cax = axCorr.scatter(self.GPS.lon, self.GPS.lat, 16, c=self.resids, cmap='viridis')
        fig.colorbar(cax, ax=axCorr, orientation='vertical')

        axCorr.set_title('Adjusted image')

        fig.tight_layout()



### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Load data
    # Instantiate object
    correction = GPScorrection(verbose=inps.verbose)

    # Load GPS data
    correction.load_gps(inps.GPSname, lonCol=inps.lonCol, latCol=inps.latCol, DATAcol=inps.DATAcol,
        headerRows=inps.headerRows, delimiter=inps.delimiter, scaleFactor=inps.scaleFactor)

    # Load image data set
    correction.load_image(inps.imgName, cohName=inps.cohName, maskArgs=inps.maskArgs)


    ## Correct to GPS
    correction.correct_to_gps(GPSclip=inps.GPSclip, fitDegree=inps.fitDegree, searchRadius=inps.searchRadius)


    ## Outputs
    # Save to file
    correction.save(inps.outName)

    # Plot if requested
    if inps.plot == True:
        correction.plot()


    plt.show()