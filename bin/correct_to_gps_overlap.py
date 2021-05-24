#!/usr/bin/env python3
'''
SHORT DESCRIPTION
Adjust multiple images to minimize the difference between a set of reference
 points (GPS) and each other, simultaneously.

FUTURE IMPROVEMENTS
    * coherence matrix
    * no sin inc norm

TESTING STATUS
In development.
'''

### IMPORT MODULES ---
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from IOsupport import append_fname, confirm_outdir, confirm_outname_ext, GPSdata, pick_dataset, load_gdal_datasets, images_from_datasets, save_gdal_dataset
from Masking import mask_datasets
from GeoFormatting import get_raster_size, lola_to_xy, xy_to_lola, transform_to_extent, transform_to_bounds
from RasterResampling import sample_points_from_raster, match_rasters
from Filtering import filter_image
from Fitting import fit_surface_to_points, design_matrix2d
from Viewing import dataset_clip_values, plot_raster


### PARSER ---
Description = '''Adjust multiple images to minimize the difference between a set of
 reference points (GPS) and each other, simultaneously.'''

Examples = ''''''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    ImgArgs = parser.add_argument_group('IMAGE INPUTS')
    ImgArgs.add_argument(dest='imgNames', type=str, nargs='+',
        help='Image file names.')
    ImgArgs.add_argument('-m','--mask', dest='maskArgs', nargs='+', type=str, default=None,
        help='Arguments for masking values/maps. ([None]).')


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
    FitArgs.add_argument('--max-overlap', dest='maxOlap', type=int, default=1000,
        help='Maximum number of overlap samples. ([1000]).')
    FitArgs.add_argument('-w','--weights', dest='weights', type=str, nargs='+', default='auto',
        help='Weights for GPS and overlap data points, respectively. Default [auto] results in GPS and overlap points being given equal status.')


    OutputArgs = parser.add_argument_group('OUTPUTS')
    OutputArgs.add_argument('-o','--out-directory', dest='outDir', type=str, default='.',
        help='Output directory.')
    OutputArgs.add_argument('-v','--verbose', dest='verbose', action='store_true',
        help='Verbose mode.')
    OutputArgs.add_argument('-p','--plot', dest='plot', action='store_true',
        help='Plot results.')

    return parser

def cmdParser(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)


### SAMPLES CLASS ---
class SampleSet:
    def __init__(self, name):
        self.name = name



### GPS CORRECTION ---
class JointCorrection:
    def __init__(self, verbose=False):
        '''
        Correct an image to GPS point measurements.
        '''
        # Parameters
        self.verbose = verbose


    ## Loading
    def load_images(self, imgNames, maskArgs=None):
        '''
        Load image data sets.
        Resample to common grid.
        Parse geographic information.
        Create mask layer.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Loading image data')

        # Load GDAL-compatible data sets
        datasets = load_gdal_datasets(imgNames, verbose=self.verbose)

        # Resample
        datasets = match_rasters(datasets, verbose=self.verbose)

        # Retrieve geographic information
        self.M, self.N = get_raster_size(pick_dataset(datasets))
        self.tnsf = pick_dataset(datasets).GetGeoTransform()
        self.extent = transform_to_extent(self.tnsf, self.M, self.N, verbose=self.verbose)
        self.bounds = transform_to_bounds(self.tnsf, self.M, self.N, verbose=self.verbose)
        self.proj = pick_dataset(datasets).GetProjection()

        self.imgs = images_from_datasets(datasets, verbose=self.verbose)

        # Create masks
        self.masks = mask_datasets(datasets, maskArgs, refBounds=self.bounds, refShape=(self.M, self.N),
            verbose=self.verbose)

        # Record parameters
        self.dsNames = list(datasets.keys())  # data set names
        self.nDatasets = len(self.dsNames)  # number of data sets


        # Report if requested
        if self.verbose == True: print('{:d} data sets loaded'.format(self.nDatasets))


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

        # Clip to map bounds
        self.GPS.crop_to_extent(*self.extent)

        # Convert GPS coordinates to pixel coordinates
        self.GPS.px, self.GPS.py = lola_to_xy(self.tnsf, self.GPS.lon, self.GPS.lat)


    ## GPS correction
    def correct_images(self, searchRadius=0, maxOlap=1000, weights='auto'):
        '''
        Correct the image to the GPS data.
        First, find the differences.
        Then, solve for the difference planes.
        Finally, add the difference planes to the images.
        '''
        # Establish grids
        self.__establish_grid__()

        # Sample GPS points
        self.__sample_gps_values__(searchRadius)

        # Sample overlapping regions
        self.__sample_overlaps__(searchRadius, maxOlap)

        # Formulate design matrix
        self.__formulate_design_matrix__()

        # Formulate weighting matrix
        self.__formulate_weighting_matrix__(weights)

        # Solve for planes
        self.__solve_for_difference_planes__()

        # Apply corrections
        self.__construct_difference_planes__()

        # Remove differences
        self.__remove_differences__()

        # Compute residuals
        self.__compute_residuals__(searchRadius)


    def __establish_grid__(self):
        '''
        Create X, Y grid for data sets.
        All data sets should be resampled to the same dimensions, so only a
         single grid is needed.
        '''
        if self.verbose == True: print('Establishing X,Y grid')

        # Establish grid
        x = np.arange(self.N)
        y = np.arange(self.M)
        self.X, self.Y = np.meshgrid(x, y)

        # Report if requested
        if self.verbose == True:
            print('grid size: {:d} x {:d}'.format(*self.X.shape))


    def __sample_gps_values__(self, searchRadius=0):
        '''
        Sample the GPS values of each image.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Sampling GPS locations using a search radius of {:d} pixels'.format(searchRadius))

        # Sample points from each image
        self.GPSsamples = {}

        for dsName in self.dsNames:
            self.GPSsamples[dsName] = SampleSet(dsName)
            self.GPSsamples[dsName].x, self.GPSsamples[dsName].y, self.GPSsamples[dsName].data, self.GPSsamples[dsName].ndx = \
                sample_points_from_raster(self.imgs[dsName], self.GPS.px, self.GPS.py, mask=self.masks[dsName], searchR=searchRadius)

            # Compute GPS - img differences
            self.GPSsamples[dsName].diff = self.GPS.data[self.GPSsamples[dsName].ndx] - self.GPSsamples[dsName].data

            # Convert pixel values back to coordinates
            self.GPSsamples[dsName].lon, self.GPSsamples[dsName].lat = xy_to_lola(self.tnsf, self.GPSsamples[dsName].x, self.GPSsamples[dsName].y)

            # Number of samples
            self.GPSsamples[dsName].n = len(self.GPSsamples[dsName].data)

            # Report
            if self.verbose == True:
                print('\t{:d} GPS samples for track {:s}'.format(self.GPSsamples[dsName].n, dsName))

        # Total number of GPS sample points
        self.nGPSsamples = np.sum([self.GPSsamples[dsName].n for dsName in self.dsNames], dtype=int)


    def __sample_overlaps__(self, searchRadius=0, maxOlap=1000):
        '''
        Sample overlaps between maps.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Sampling overlapping regions of images')

        # Average images
        if self.verbose == True:
            print('Resampling images with width {:d} boxcar'.format(searchRadius))

        w = int(2*searchRadius + 1)  # window width
        avedImgs = [filter_image('mean', self.imgs[dsName], w) for dsName in self.dsNames]

        # Find areas of overlap
        self.OLAPnames = []
        self.OLAPsamples = {}

        for i in range(self.nDatasets-1):
            dsName1 = self.dsNames[i]

            for j in range(i+1, self.nDatasets):
                dsName2 = self.dsNames[j]

                # Check whether images overlap
                olap = (self.masks[dsName1]==1) & (self.masks[dsName2]==1)
                nOlap = np.sum(olap)  # number of overlapping points

                # Proceed if there is overlap between images
                if nOlap > 0:
                    # Formulate name of overlap region
                    OLAPname = '{:s}_{:s}'.format(dsName1, dsName2)
                    self.OLAPnames.append(OLAPname)

                    # Instantiate sample set
                    self.OLAPsamples[OLAPname] = SampleSet(OLAPname)

                    # Indices of overlapping images
                    self.OLAPsamples[OLAPname].ndx1 = i
                    self.OLAPsamples[OLAPname].ndx2 = j

                    # Report if requested
                    if self.verbose == True:
                        print('{:d} overlap points between tracks {:s} (ndx {:d}) and {:s} (ndx {:d})'.\
                            format(nOlap, dsName1, i, dsName2, j))

                    # Determine downsample factor if necessary
                    if nOlap > maxOlap:
                        skips = int(np.ceil(nOlap/maxOlap))
                        if self.verbose == True: print('enforcing max samples {:d}'.format(maxOlap))
                    else:
                        skips = int(1)

                    # Sample pixel coordinates
                    x = self.X[olap].flatten()
                    y = self.Y[olap].flatten()

                    self.OLAPsamples[OLAPname].x = x[::skips]
                    self.OLAPsamples[OLAPname].y = y[::skips]

                    # Sample differences
                    diff = (avedImgs[j][olap] - avedImgs[i][olap]).flatten()

                    self.OLAPsamples[OLAPname].diff = diff[::skips]

                    # Number of samples
                    self.OLAPsamples[OLAPname].n = len(self.OLAPsamples[OLAPname].diff)

        # Total number of overlap samples
        self.nOLAPsamples = np.sum([self.OLAPsamples[OLAPname].n for OLAPname in self.OLAPnames], dtype=int)


    def __formulate_design_matrix__(self):
        '''
        Formulate the design matrix including all GPS and overlappinig regions.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Constructing design matrix')

        # Total number of sample points
        self.nTotalSamples = self.nGPSsamples + self.nOLAPsamples
        
        # Report if requested
        if self.verbose == True:
            print('{:d} total GPS samples'.format(self.nGPSsamples))
            print('{:d} total OLAP samples'.format(self.nOLAPsamples))
            print('{:d} total samples'.format(self.nTotalSamples))

        # Construct design matrix
        self.G = np.zeros((self.nTotalSamples, 3*self.nDatasets))

        # Construct observation vector
        self.diffs = np.zeros((self.nTotalSamples, 1))

        # Populate design matrix and observations with GPS samples
        startNdx = 0
        for i, dsName in enumerate(self.dsNames):
            # Number of sample points in this image
            nSamples = self.GPSsamples[dsName].n

            # Fill in design matrix
            self.G[startNdx:startNdx+nSamples, i*3+0] = 1
            self.G[startNdx:startNdx+nSamples, i*3+1] = self.GPSsamples[dsName].x
            self.G[startNdx:startNdx+nSamples, i*3+2] = self.GPSsamples[dsName].y

            # Fill in differences vector
            self.diffs[startNdx:startNdx+nSamples, 0] = self.GPSsamples[dsName].diff

            # Update counter
            startNdx += nSamples

        # Populate design matrix and observatoins with overlap samples
        for OLAPname in self.OLAPnames:
            # Number of samples points in this pair of images
            nSamples = self.OLAPsamples[OLAPname].n

            # Matrix indices
            ndx1 = 3*self.OLAPsamples[OLAPname].ndx1
            ndx2 = 3*self.OLAPsamples[OLAPname].ndx2

            # Fill in the design matrix
            self.G[startNdx:startNdx+nSamples, ndx1+0] = 1
            self.G[startNdx:startNdx+nSamples, ndx1+1] = self.OLAPsamples[OLAPname].x
            self.G[startNdx:startNdx+nSamples, ndx1+2] = self.OLAPsamples[OLAPname].y
            self.G[startNdx:startNdx+nSamples, ndx2+0] = -1
            self.G[startNdx:startNdx+nSamples, ndx2+1] = -self.OLAPsamples[OLAPname].x
            self.G[startNdx:startNdx+nSamples, ndx2+2] = -self.OLAPsamples[OLAPname].y

            # Fill in the differences vector
            self.diffs[startNdx:startNdx+nSamples, 0] = self.OLAPsamples[OLAPname].diff

            # Update counter
            startNdx += nSamples


    def __formulate_weighting_matrix__(self, weights='auto'):
        '''
        Formulate the identity-like matrix of relative GPS and overlap weights.
        Setting the weighting to 'auto' means the GPS overlap samples will be
         weighted the same in the final scheme.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Constructing weighting matrix')

        # Determine weights
        if weights == 'auto':
            GPSweight = np.ceil(self.nTotalSamples/self.nGPSsamples)
            OLAPweight = np.ceil(self.nTotalSamples/self.nOLAPsamples)

            if self.verbose == True:
                print('Automatically adjusting weights to GPS: {:f}; overlap: {:f}'.format(GPSweight, OLAPweight))

        elif len(weights) == 2:
            GPSweight = float(weights[0])
            OLAPweight = float(weights[1])

            if self.verbose == True:
                print('Using GPS weight: {:f}; overlap weight: {:f}'.format(GPSweight, OLAPweight))

        else:
            print('Cannot identify weighting values')
            exit()

        # Weighting matrix
        weights = np.ones(self.nTotalSamples)
        weights[:self.nGPSsamples] = GPSweight
        weights[self.nGPSsamples:self.nOLAPsamples] = OLAPweight
        self.W = np.identity(self.nTotalSamples)
        for i in range(self.nTotalSamples):
            self.W[i,i] = weights[i]


    def __solve_for_difference_planes__(self):
        '''
        Invert design matrix and use data to solve for best-fit planes to the
         differences.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Solving for planes')

        # Invert design matrix
        Ginv = np.linalg.inv(self.G.T.dot(self.W).dot(self.G)).dot(self.G.T).dot(self.W)

        # Solve for planes
        self.B = Ginv.dot(self.diffs)


    def __construct_difference_planes__(self):
        '''
        Compute the difference maps from the inversion solution.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Computing difference planes')

        # Design matrix for difference planes
        Gdiff = design_matrix2d(self.X.flatten(), self.Y.flatten(), degree=1)

        # Construct difference planes
        self.diffPlanes = {}
        for i, dsName in enumerate(self.dsNames):
            i3 = i*3
            self.diffPlanes[dsName] = Gdiff.dot(self.B[i3:i3+3, 0]).reshape(self.M, self.N)


    def __remove_differences__(self):
        '''
        Add the difference plane to each image.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Adding difference planes to original data')

        # Correct images
        self.corrImgs = {}
        for dsName in self.dsNames:
            self.corrImgs[dsName] = self.imgs[dsName] + self.diffPlanes[dsName]


    def __compute_residuals__(self, searchRadius=0):
        '''
        Compute data - GPS residuals.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Computing residuals')

        # Data - GPS residuals
        self.GPSresids = {}
        for dsName in self.dsNames:
            self.GPSresids[dsName] = SampleSet(dsName)

            # Sample corrected data points
            self.GPSresids[dsName].x, self.GPSresids[dsName].y, self.GPSresids[dsName].corrData, self.GPSresids[dsName].ndx = \
                sample_points_from_raster(self.corrImgs[dsName], self.GPS.px, self.GPS.py,
                    mask=self.masks[dsName], searchR=searchRadius, verbose=False)

            # Convert pixel values to lon lat coordinates
            self.GPSresids[dsName].lon, self.GPSresids[dsName].lat = xy_to_lola(self.tnsf, self.GPSresids[dsName].x, self.GPSresids[dsName].y)

            # Compute residuals
            self.GPSresids[dsName].resids = self.GPSresids[dsName].corrData - self.GPS.data[self.GPSresids[dsName].ndx]


    ## Saving
    def save(self, outDir):
        '''
        Save image to file.
        '''
        if inps.verbose == True:
            print('*'*32)
            print('Saving outputs')

        # Confirm that output directory exists
        outDir = os.path.abspath(outDir)

        if not os.path.exists(outDir):
            os.mkdir(outDir)

        # Loop through images
        for dsName in self.dsNames:
            outName = os.path.join(outDir, '{:s}_corrected.tif'.format(dsName))

            # Save adjusted images to file
            save_gdal_dataset(outName, self.corrImgs[dsName], mask=self.masks[dsName],
                proj=self.proj, tnsf=self.tnsf,
                verbose=self.verbose)


    ## Plotting
    def plot(self):
        '''
        Plot inputs and outputs.
        '''
        if self.verbose == True:
            print('*'*32)
            print('Plotting results')

        # Plot inputs
        self.__plot_inputs__()

        # Plot differences
        self.__plot_differences__()

        # Plot results
        self.__plot_corrected_maps__()

        # Velocity vs GPS plot
        self.__plot_img_vs_gps__()


    def __plot_inputs__(self):
        '''
        Plot input data sets.
        '''
        # Spawn figure and axis
        OGfig, axOG = plt.subplots(figsize=(8,7))

        # Color parameters
        minPct = 1
        maxPct = 99
        cbarOrient = 'auto'

        # Colorscale
        mapsMin, mapsMax = dataset_clip_values(self.imgs, minPct=minPct, maxPct=maxPct, masks=list(self.masks.values()))
        gpsMin, gpsMax = np.percentile(self.GPS.data, [minPct, maxPct])
        vmin = np.min([mapsMin, gpsMin])
        vmax = np.max([mapsMax, gpsMax])

        # Plot rasters
        for i, dsName in enumerate(self.dsNames):
            if i == 0:
                cbarOrient = cbarOrient
            else:
                cbarOrient = None

            plot_raster(self.imgs[dsName], extent=self.extent, mask=self.masks[dsName],
                vmin=vmin, vmax=vmax, cbarOrient=cbarOrient,
                fig=OGfig, ax=axOG)

        # Plot GPS
        axOG.scatter(self.GPS.lon, self.GPS.lat, 16, c=self.GPS.data,
            cmap='viridis', vmin=vmin, vmax=vmax)

        # Format plot
        axOG.set_title('Input data')


    def __plot_differences__(self):
        '''
        Plot difference ramps and GPS differences.
        '''
        # Spawn figure and axis
        DiffFig, axDiff = plt.subplots(figsize=(8,7))

        # Color parameters
        minPct = 1
        maxPct = 99
        cbarOrient = 'auto'

        # Colorscale
        gpsDiffMins = []; gpsDiffMaxs = []
        for dsName in self.dsNames:
            gpsDiffMins.append(self.GPSsamples[dsName].diff.min())
            gpsDiffMaxs.append(self.GPSsamples[dsName].diff.max())

        vmin, vmax = np.min(gpsDiffMins), np.max(gpsDiffMaxs)

        # Plot difference planes
        for i, dsName in enumerate(self.dsNames):
            if i == 0:
                cbarOrient = cbarOrient
            else:
                cbarOrient = None

            plot_raster(self.diffPlanes[dsName], extent=self.extent, mask=self.masks[dsName],
                cmap='jet', vmin=vmin, vmax=vmax, cbarOrient=cbarOrient,
                fig=DiffFig, ax=axDiff)

        # Plot GPS differences
        for dsName in self.dsNames:
            axDiff.scatter(self.GPSsamples[dsName].lon, self.GPSsamples[dsName].lat, 16, c=self.GPSsamples[dsName].diff,
                cmap='jet', vmin=vmin, vmax=vmax)

        # Format plot
        axDiff.set_title('Differences')
        axDiff.set_aspect(1)


    def __plot_corrected_maps__(self):
        '''
        Plot adjusted images and residuals.
        '''
        # Spawn figure and axis
        CorrFig, axCorr = plt.subplots(figsize=(8,8))

        # Color parameters
        minPct = 1
        maxPct = 99
        cbarOrient = 'horizontal'

        # Colorscales
        mapsMin, mapsMax = dataset_clip_values(self.corrImgs, minPct=minPct, maxPct=maxPct, masks=list(self.masks.values()))
        gpsMins, gpsMaxs = [], []
        for dsName in self.dsNames:
            gpsMin, gpsMax = np.percentile(self.GPSresids[dsName].resids, (minPct, maxPct))
            gpsMins.append(gpsMin)
            gpsMaxs.append(gpsMax)
        gpsMin = np.min(gpsMins)
        gpsMax = np.max(gpsMaxs)

        # Plot rasters
        for i, dsName in enumerate(self.dsNames):
            if i == 0:
                cbarOrient = cbarOrient
            else:
                cbarOrient = None

            plot_raster(self.corrImgs[dsName], extent=self.extent, mask=self.masks[dsName],
                vmin=mapsMin, vmax=mapsMax, cbarOrient=cbarOrient,
                fig=CorrFig, ax=axCorr)

        # Plot GPS residuals
        for i, dsName in enumerate(self.dsNames):
            cax = axCorr.scatter(self.GPSresids[dsName].lon, self.GPSresids[dsName].lat, 16, c=self.GPSresids[dsName].resids,
                cmap='jet', vmin=gpsMin, vmax=gpsMax)

            if i == 0:
                CorrFig.colorbar(cax, ax=axCorr, orientation='vertical')


    def __plot_img_vs_gps__(self):
        '''
        Plot the image values at each GPS location vs the GPS value at that
         location.

        Corrected image values are stored in the GPSresids objects.
        '''
        # Spawn figure and axis
        GPSfig, axGPS = plt.subplots()

        # Reference values
        GPSmin, GPSmax = [], []

        # Plot corrected image values vs GPS
        for dsName in self.dsNames:
            # Retrieve values
            gps = self.GPS.data[self.GPSresids[dsName].ndx]
            data = self.GPSresids[dsName].corrData

            # Plot GPS and data
            axGPS.scatter(gps, data, label=dsName)

            # Append min/max values
            GPSmin.append(np.min(gps))
            GPSmax.append(np.max(gps))

        axGPS.legend()

        # Plot 1:1 line
        GPSmin = np.min(GPSmin)*1.2
        GPSmax = np.max(GPSmax)*1.2

        axGPS.plot([GPSmin, GPSmax], [GPSmin, GPSmax], 'k--')




### MAIN ---
if __name__ == '__main__':
    ## Inputs
    # Gather arguments
    inps = cmdParser()


    ## Load data
    # Instantiate object
    correction = JointCorrection(verbose=inps.verbose)

    # Load image data set
    correction.load_images(inps.imgNames, maskArgs=inps.maskArgs)

    # Load GPS data
    correction.load_gps(inps.GPSname, lonCol=inps.lonCol, latCol=inps.latCol, DATAcol=inps.DATAcol,
        headerRows=inps.headerRows, delimiter=inps.delimiter, scaleFactor=inps.scaleFactor)


    ## Correct to GPS and minimize overlap
    correction.correct_images(inps.searchRadius, inps.maxOlap, inps.weights)


    ## Outputs
    # Save to file
    correction.save(inps.outDir)

    # Plot if requested
    if inps.plot == True:
        correction.plot()


    plt.show()