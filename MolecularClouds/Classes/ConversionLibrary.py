import numpy as np


def ra_hms2deg(ra_h, ra_m, ra_s):
    """
     This function converts a right ascension in hour:min:sec to degrees
    :param ra_h: hour component of right ascension
    :param ra_m: minute component of right ascension
    :param ra_s: second component of right ascension
    :return: Right ascension in degrees
    """
    return np.array(ra_h) * 15 + np.array(ra_m)/4 + np.array(ra_s)/240


def dec_dms2deg(dec_d, dec_m, dec_s):
    """
     This function converts a declination in degree:arcmin:arcsec to degrees
    :param dec_d: degree component of declination
    :param dec_m: arcminute component of declination
    :param dec_s: arcsecond component of declination
    :return: Declination in degrees
    """
    return (np.array(abs(dec_d)) + np.array(abs(dec_m))/60 + np.array(abs(dec_s))/3600) * np.sign(dec_d)


def ra_deg2hms(ra_deg):
    """
     This function converts a right ascension in degrees to hour:min:sec
    :param: ra_deg: Right ascension in degrees
    :return ra_h: hour component of right ascension
    :return ra_m: minute component of right ascension
    :return ra_s: second component of right ascension
    """
    ra_h = ra_deg // 15
    remaining_deg = ra_deg % 15
    ra_m = remaining_deg // (1/4)
    remaining_deg = remaining_deg % (1/4)
    ra_s = remaining_deg // (1/240)

    return ra_h, ra_m, ra_s


def dec_deg2dms(dec_deg):
    """
     This function converts a declination in degree:arcmin:arcsec to degrees

    :param: dec_deg: Declination in degrees
    :return dec_d: degree component of declination
    :return dec_m: arcminute component of declination
    :return dec_s: arcsecond component of declination
    """

    dec_d = dec_deg // 1
    remaining_deg = dec_deg % 1
    dec_m = remaining_deg // (1/60)
    remaining_deg = remaining_deg % (1/60)
    dec_s = remaining_deg // (1/3600)

    return dec_d, dec_m, dec_s


def RADec2xy(RA, Dec, wcs):
    xCoords = []
    yCoords = []
    for i in range(len(RA)):
        pixelRow, pixelColumn = wcs.wcs_world2pix(RA[i], Dec[i], 0)
        xCoords.append(pixelRow)
        yCoords.append(pixelColumn)
    return xCoords, yCoords
