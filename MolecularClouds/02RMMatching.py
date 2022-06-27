"""
This is the second stage of the BLOSMapping method where rotation measure points from the catalogue are matched to
an extinction value from the fits file. Matching is based on physical proximity.

The matched rotation measure data and extinction information are saved in a file.
"""
import os
import csv

from itertools import zip_longest

from astropy.io import fits
from astropy.wcs import WCS

import numpy as np
import math

from Classes.DataFile import DataFile
from Classes.RegionOfInterest import Region
import Classes.config as config

from Classes.util import getBoxBounds
import Classes.ConversionLibrary as cl
import Classes.RefJudgeLib as rjl

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = input("Enter the name of the region of interest: ")
cloudName = cloudName.capitalize()  # Ensure only the first letter is capitalized
#cloudName = config.cloud
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
RMCatalogPath = os.path.join(config.dir_root, config.dir_data, config.file_RMCatalogue)
saveFilePath = os.path.join(config.dir_root, config.dir_fileOutput, cloudName, config.prefix_RMExtinctionMatch + cloudName + '.txt')
# -------- DEFINE FILES AND PATHS. --------

# -------- READ FITS FILE --------
hdulist = fits.open(regionOfInterest.fitsFilePath)
hdu = hdulist[0]
wcs = WCS(hdu.header)
# -------- READ FITS FILE. --------

# -------- PREPROCESS FITS DATA TYPE. --------
data = rjl.deepCopy(hdu.data)

# If fitsDataType is column density, then convert to visual extinction
if regionOfInterest.fitsDataType == 'HydrogenColumnDensity':
    data = data / config.VExtinct_2_Hcol

#Handle bad data (negative/no values) by interpolation.
boxXMin = regionOfInterest.xmin
boxXMax = regionOfInterest.xmax
boxYMin = regionOfInterest.ymin
boxYMax = regionOfInterest.ymax

xmin, xmax, ymin, ymax = getBoxBounds(data, boxXMin, boxXMax, boxYMin, boxYMax)

data[ymin:ymax, xmin:xmax][data[ymin:ymax, xmin:xmax] < 0] = np.nan
baddata = np.isnan(data)
if config.doInterpExtinct:
    data[ymin:ymax, xmin:xmax] = rjl.interpMask(data[ymin:ymax, xmin:xmax], baddata[ymin:ymax, xmin:xmax], 'linear') #This step is computationally costly. It may be omitted if it is taking too long.
# -------- PREPROCESS FITS DATA TYPE. --------

# -------- READ ROTATION MEASURE FILE --------
# Get all the rm points within the region of interest
rmData = DataFile(RMCatalogPath, regionOfInterest.raHoursMax, regionOfInterest.raMinsMax, regionOfInterest.raSecMax,
                   regionOfInterest.raHoursMin, regionOfInterest.raMinsMin, regionOfInterest.raSecMin,
                   regionOfInterest.decDegMax, regionOfInterest.decDegMin)
# -------- READ ROTATION MEASURE FILE. --------

# -------- DEFINE THE ERROR RANGE --------
# The physical limit on how far an extinction value can be from the rm and still be considered valid/applicable
# Uncertainty based.
raErrsSec = np.array(rmData.targetRAErrSecs)
decErrs = np.array(rmData.targetDecErrArcSecs)

raErrSec = max(abs(raErrsSec)) #s
raErr = cl.ra_hms2deg(0, 0, raErrSec) #deg
decErrSec = max(abs(decErrs)) #s
decErr = cl.dec_dms2deg(0, 0, decErrSec) #deg

RMResolutionDegs = max(raErr, decErr)

ExtinctionResolutionDegs = min(abs(hdu.header['CDELT1']), abs(hdu.header['CDELT2'])) #deg
# -------- It is 1 pixel at most if the extinction map has a lower resolution than the RM map. The maximum number of pixels which fit within the RM's resolution otherwise.
if (ExtinctionResolutionDegs > RMResolutionDegs):
    NDelt = 1
else:
    NDelt = np.ceil(RMResolutionDegs/ExtinctionResolutionDegs)
# -------- DEFINE THE ERROR RANGE. --------

# -------- DEFINE PARAMETERS --------
ExtinctionIndex_x = []
ExtinctionIndex_y = []
ExtinctionRa = []
ExtinctionDec = []
ExtinctionValue = []

Identifier = []
RMRa = []
RMDec = []
RMValue = []
RMErr = []

ErrRangePix = []

Extinction_MinInRangeRa = []
Extinction_MinInRangeDec = []
Extinction_MinInRange = []

Extinction_MaxInRangeRa = []
Extinction_MaxInRangeDec = []
Extinction_MaxInRange = []

IsExtinctionObserved = []
# -------- DEFINE PARAMETERS. --------

# -------- MATCH ROTATION MEASURES AND EXTINCTION VALUES --------
cntr = 0  # To keep track of how many matches have been made - numbering starts at 0
# Go through all of the rotation measure values and match them to an extinction value
for index in range(len(rmData.targetRotationMeasures)):

    # ---- Location of the rotation measure
    rmRA = rmData.targetRaHourMinSecToDeg[index]
    rmDec = rmData.targetDecDegArcMinSecs[index]
    py, px = wcs.world_to_array_index_values(rmRA, rmDec)  # Array indices of the rotation measure
    # ---- Location of the rotation measure.

    # If the rm lies within the given fits file:
    if 0 <= px < data.shape[1] and 0 <= py < data.shape[0]:

        # If the rm lies on a point with data:
        if data[py, px] != -1 and math.isnan(data[py, px]) is False:
            extinction = data[py, px]

            Identifier.append(cntr)
            RMRa.append(rmRA)
            RMDec.append(rmDec)
            RMValue.append(rmData.targetRotationMeasures[index])
            RMErr.append(rmData.targetRMErrs[index])

            # ---- Match rotation measure to an extinction value
            ExtinctionIndex_x.append(int(px))
            ExtinctionIndex_y.append(int(py))
            ExtinctionRa.append(wcs.wcs_pix2world(px, py, 0)[0])
            ExtinctionDec.append(wcs.wcs_pix2world(px, py, 0)[1])
            ExtinctionValue.append(extinction)
            # ---- Match rotation measure to an extinction value.

            # ---- Find the extinction error range for the given rm
            ErrRangePix.append(NDelt)
            ind_xmin, ind_xmax, ind_ymin, ind_ymax = rjl.getBoxRange(px, py, data, NDelt)
            # ---- Find the extinction error range for the given rm.

            # ---- Cycle through extinction values within the error range
            extinction_temp = []
            ra_temp = []
            dec_temp = []

            for pxx in range(ind_xmin, ind_xmax):
                for pyy in range(ind_ymin, ind_ymax):
                    extinction = data[pyy, pxx]
                    extinction_temp.append(extinction)
                    xx, yy = wcs.wcs_pix2world(pxx, pyy, 0)
                    ra_temp.append(xx)
                    dec_temp.append(yy)
            # ---- Cycle through extinction values within the error range.

            # Find minimum extinction value
            ind_min = np.where(extinction_temp == min(extinction_temp))[0][0]
            Extinction_MinInRangeRa.append(ra_temp[ind_min])
            Extinction_MinInRangeDec.append(dec_temp[ind_min])
            Extinction_MinInRange.append(extinction_temp[ind_min])

            # Find maximum extinction value
            ind_max = np.where(extinction_temp == max(extinction_temp))[0][0]
            Extinction_MaxInRangeRa.append(ra_temp[ind_max])
            Extinction_MaxInRangeDec.append(dec_temp[ind_max])
            Extinction_MaxInRange.append(extinction_temp[ind_max])
            cntr += 1

            # ---- Negative extinction (the rm value landed on a negative pixel)
            # Negative extinction is not physical; in prior step it was interpolated away. Mark these points.
            if baddata[py, px]:
                IsExtinctionObserved.append(False)
            else:
                IsExtinctionObserved.append(True)
            # ---- Negative extinction.
# -------- MATCH ROTATION MEASURES AND EXTINCTION VALUES. --------

# -------- WRITE TO A FILE --------
with open(saveFilePath, 'w') as f:
    f.write('ID#' + '\t' + 'Extinction_Index_x' + '\t' + 'Extinction_Index_y' + '\t' + 'Ra(deg)' + '\t' +
            'Dec(deg)' + '\t' + 'Rotation_Measure(rad/m2)' + '\t' + 'RM_Err(rad/m2)' + '\t' +
            'RA_inExtincFile(degree)' + '\t' +
            'Dec_inExtincFile(degree)' + '\t' + 'Extinction_Value' + '\t' + 'Error_Range(pix)' + '\t' +
            'Min_Extinction_Value' + '\t' + 'Min_Extinction_Ra' + '\t' + 'Min_Extinction_Dec' + '\t' +
            'Max_Extinction_Value' + '\t' + 'Max_Extinction_RA' + '\t' + 'Max_Extinction_dec' + '\n')
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip_longest(Identifier, ExtinctionIndex_x, ExtinctionIndex_y, RMRa, RMDec, RMValue, RMErr,
                                 ExtinctionRa, ExtinctionDec, ExtinctionValue, ErrRangePix,
                                 Extinction_MinInRange, Extinction_MinInRangeRa, Extinction_MinInRangeDec,
                                 Extinction_MaxInRange, Extinction_MaxInRangeRa, Extinction_MaxInRangeDec,
                                 fillvalue=''))
# -------- WRITE TO A FILE. --------

print('\nWithin the specified region of interest, a total of {} rotation measure points were matched '
      'to visual extinction values.\n'.format(len(Identifier)))
print('Matched visual extinction and rotation measure data were saved to {}'.format(saveFilePath))
