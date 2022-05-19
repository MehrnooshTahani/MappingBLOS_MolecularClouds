"""
This file contains the class definition needed to read rotation measure data from the Taylor et al (2009)
catalogue

-  Use of the SkyCoord package to convert coordinates may increase runtime. Conversions may be attempted with the
    SkyCoord package if desired.  Code to accomplish this is left in comments throughout the file
"""
import numpy as np
from astropy.coordinates import SkyCoord


# -------- FUNCTION DEFINITION --------
def ra_hms2deg(ra_h, ra_m, ra_s):
    """
     This function converts a right ascension in hour:min:sec to degrees
    :param ra_h: hour component of right ascension
    :param ra_m: minute component of right ascension
    :param ra_s: second component of right ascension
    :return: Right ascension in degrees
    """
    if (ra_h < 0 or ra_h > 24) or (ra_m < 0 or ra_m >= 60) or (ra_s < 0 or ra_s >= 60):
        print('Invalid RA')
        return None
    else:
        ra_deg = (ra_h * 15) + (ra_m / 4) + (ra_s / 240)
    return ra_deg
# -------- FUNCTION DEFINITION. --------


# -------- FUNCTION DEFINITION --------
def dec_dms2deg(dec_d, dec_m, dec_s):
    """
     This function converts a declination in degree:arcmin:arcsec to degrees
    :param dec_d: degree component of declination
    :param dec_m: arcminute component of declination
    :param dec_s: arcsecond component of declination
    :return: Declination in degrees
    """
    if dec_d < 0:
        decdeg = -(abs(dec_d) + (abs(dec_m)/60) + (abs(dec_s)/3600))
    else:
        decdeg = abs(dec_d) + (abs(dec_m)/60) + (abs(dec_s)/3600)

    return decdeg
# -------- FUNCTION DEFINITION. --------


# -------- CLASS DEFINITION --------
class DataFile:
    def __init__(self, filename, raHoursMax, raMinsMax, raSecMax, raHoursMin, raMinsMin,
                 raSecMin, decDegMax, decDegMin):
        """
           Takes a file containing rotation measure data in the format of the Taylor et al. (2009) catalog and gives
           parameters such as ra, dec, rm, etc corresponding to a specific region of interest

        :param filename: Path to the file containing the rotation measure data (eg the Taylor et al (2009) catalogue)
        :param raHoursMax: Hour component of the maximum right ascension of the region of interest
        :param raMinsMax:  Minute component of the maximum right ascension of the region of interest
        :param raSecMax:   Second component of the maximum right ascension of the region of interest
        :param raHoursMin: Hour component of the minimum right ascension of the region of interest
        :param raMinsMin:  Minute component of the minimum right ascension of the region of interest
        :param raSecMin:   Second component of the minimum right ascension of the region of interest
        :param decDegMax:  Degrees of the maximum declination of the region of interest
        :param decDegMin:  Degrees of the minimum declination of the region of interest

        Notes
        --------

        Structure of the Taylor et al. (2009) Rotation Measure Catalogue
        --------------------------------------------------------------------------------
        columns  Format  Units   Quantity
        --------------------------------------------------------------------------------
          1-  2 I2      h       Hour of Right Ascension (J2000)
          4-  5 I2      min     Minute of Right Ascension (J2000)
          7- 11 F5.2    s       Second of Right Ascension (J2000)
         13- 16 F4.2    s       1-sigma error in Right Ascension
         18- 20 I3      deg     Degree of Declination (J2000)
         22- 23 I2      arcmin  Arcminute of Declination (J2000)
        25- 29 F5.2    arcsec  Arcsecond of Declination (J2000)
        31- 33 F3.1    arcsec  1-sigma error in Declination
        35- 43 F9.4    deg     Galactic longitude
        45- 52 F8.4    deg     Galactic latitude
        54- 60 F7.1    mJy     NVSS integrated Stokes I flux density
        62- 66 F5.1    mJy     1-sigma error in Stokes I flux density
        68- 74 F7.2    mJy     Average peak polarized intensity
        76- 80 F5.2    mJy     1-sigma error in polarized intensity
        82- 87 F6.2    %       Percent polarization (1)
        89- 93 F5.2    %       1-sigma error in m
        95-100 F6.1   rad/m2   Rotation Measusure
        102-105 F4.1   rad/m2   1-sigma error in RM
        """

        raHours, raMins, raSecs, raErrSecs, decDegs, decArcmins, decArcsecs, decErrArcsecs, longitudeDegs, \
            latitudeDegs, nvssStokesIs, stokesIErrs, AvePeakPIs, PIErrs, polarizationPercets, mErrPercents, \
            rotationMeasures, RMErrs = np.loadtxt(filename, usecols=(0, 1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 14, 15, 17,
                                                                     18, 20, 21, 23), unpack=True)

        self.targetRAHours = []  # Hour component of right ascension in hr:min:sec
        self.targetRAMins = []  # Minute component of right ascension in hr:min:sec
        self.targetRASecs = []  # Second component of right ascension in hr:min:sec
        self.targetRaMinsSecs = []
        self.targetRaHourMinSecToDeg = []
        self.targetRAErrSecs = []
        self.targetDecDegs = []  # Degree component of declination in deg:arcmin:arcsec
        self.targetDecArcMins = []  # Arcminute component of declination in deg:arcmin:arcsec
        self.targetDecArcSecs = []  # Arcsecond component of declination in deg:arcmin:arcsec
        self.targetDecErrArcSecs = []
        self.targetDecDegArcMinSecs = []
        self.targetLongitudeDegs = []
        self.targetLatitudeDegs = []
        self.targetNvssStokesIs = []
        self.targetStokesIErrs = []
        self.targetAvePeakPIs = []
        self.targetPIErrs = []
        self.targetPolarizationPercets = []
        self.targetMErrPercents = []
        self.targetRotationMeasures = []
        self.targetRMErrs = []

        radegMax = ra_hms2deg(raHoursMax, raMinsMax, raSecMax)  # <- If converting manually
        radegMin = ra_hms2deg(raHoursMin, raMinsMin, raSecMin)  # <- If converting manually

        # If using SkyCoord:
        # -------- MAKE SKYCOORD OBJECTS FOR REGION OF INTEREST  --------
        # ra_hmsMax = str(int(raHoursMax))+'h'+str(int(raMinsMax))+'m'+str(raSecMax)+'s'
        # coord_Max = SkyCoord(ra=ra_hmsMax, dec=decDegMax, unit='degree')
        # radegMax = coord_Max.ra.degree
        #
        # ra_hmsMin = str(int(raHoursMin))+'h'+str(int(raMinsMin))+'m'+str(raSecMin)+'s'
        # coord_Min = SkyCoord(ra=ra_hmsMin, dec=decDegMin, unit='degree')
        # radegMin = coord_Min.ra.degree
        # -------- MAKE SKYCOORD OBJECTS FOR REGION OF INTEREST.  --------

        # -------- EXTRACT INFORMATION FROM THE REGION OF INTEREST --------
        for index in range(len(raHours)):

            # If converting manually:
            radeg = ra_hms2deg(raHours[index], raMins[index], raMins[index])  # <- If converting manually
            decdeg = dec_dms2deg(decDegs[index], decArcmins[index], decArcsecs[index])  # <- If converting manually

            # If using SkyCoord:
            # ---- Make given coordinate into a SkyCoord object and convert to degrees
            # ra_hms = str(int(raHours[index])) + 'h' + str(int(raMins[index])) + 'm' + str(raSecs[index]) + 's'
            # dec_dms = str(int(decDegs[index])) + 'd' + str(int(decArcmins[index])) + 'm' + str(decArcsecs[index]) + 's'
            # coord = SkyCoord(ra=ra_hms, dec=dec_dms)
            # radeg = coord.ra.degree
            # decdeg = coord.dec.degree
            # ---- Make given coordinate into a SkyCoord object and convert to degrees.

            # If the coordinate is within the region of interest, extract information
            if radegMin <= radeg < radegMax and decDegMin <= decdeg <= decDegMax:
                self.targetRAHours.append(raHours[index])
                self.targetRAMins.append(raMins[index])
                self.targetRaMinsSecs.append(raMins[index] + ((raSecs[index]) / 60.0))
                self.targetRaHourMinSecToDeg.append(radeg)
                self.targetRASecs.append(raSecs[index])
                self.targetRAErrSecs.append(raErrSecs[index])
                self.targetDecDegs.append(decDegs[index])
                self.targetDecArcMins.append(decArcmins[index])
                self.targetDecArcSecs.append(decArcsecs[index])
                self.targetDecDegArcMinSecs.append(decdeg)
                self.targetDecErrArcSecs.append(decErrArcsecs[index])
                self.targetLongitudeDegs.append(longitudeDegs[index])
                self.targetLatitudeDegs.append(latitudeDegs[index])
                self.targetNvssStokesIs.append(nvssStokesIs[index])
                self.targetStokesIErrs.append(stokesIErrs[index])
                self.targetAvePeakPIs.append(AvePeakPIs[index])
                self.targetPIErrs.append(PIErrs[index])
                self.targetPolarizationPercets.append(polarizationPercets[index])
                self.targetMErrPercents.append(mErrPercents[index])
                self.targetRotationMeasures.append(rotationMeasures[index])
                self.targetRMErrs.append(RMErrs[index])
        # -------- EXTRACT INFORMATION FROM THE REGION OF INTEREST. --------
# -------- CLASS DEFINITION. --------
