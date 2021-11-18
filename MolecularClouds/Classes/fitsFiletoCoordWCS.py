from astropy.io import fits
from astropy.wcs import WCS

"""
This file contains a class which reads a fits file and converts the pixels to equatorial coordinate values 
and makes an ndarray.
"""


# -------- CLASS DEFINITION --------
class FitsFileToCoordWCS:
    def __init__(self, fitsfilepath):
        hdulist = fits.open(fitsfilepath)
        hdu = hdulist[0]
        wcs = WCS(hdu.header)

        self.densityRA = []
        self.densityDec = []
        self.extinctionValue = []

        for i in range(hdu.data.shape[1]):
            for j in range(hdu.data.shape[0]):
                px, py = wcs.wcs_pix2world(i, j, 0)
                value = hdu.data[j, i]
                if value >= 0:
                    self.densityRA.append(px)
                    self.densityDec.append(py)
                    self.extinctionValue.append(value)

        hdulist.close()
# -------- CLASS DEFINITION. --------
