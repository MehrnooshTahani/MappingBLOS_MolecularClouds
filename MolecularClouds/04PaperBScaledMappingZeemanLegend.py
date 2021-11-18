"""
This file is the fourth step of the BLOSMapping method where the line of sight component of the magnetic field is
plotted as scatter points on the extinction map.  This file will also plot Zeeman measurements as scatter points.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
import os
from MolecularClouds.Classes.RegionOfInterest import Region

# -------- CHOOSE THE REGION OF INTEREST --------
cloudName = 'California'
# cloudName = input("Enter the name of the region of interest: ")
cloudName = cloudName.capitalize()  # Ensure only the first letter is capitalized
regionOfInterest = Region(cloudName)
# -------- CHOOSE THE REGION OF INTEREST. --------

# -------- DEFINE FILES AND PATHS --------
currentDir = os.path.abspath(os.getcwd())
BScaledFilePath = os.path.join(currentDir, 'FileOutput/BScaledNRefCalifornia12.txt')
saveFilePath = os.path.join(currentDir, 'ResultPlotPDF/BMapping'+cloudName+'.png')
# -------- DEFINE FILES AND PATHS. --------

# -------- LOAD BLOS DATA --------
raDegree, decDegree, BScaled = np.loadtxt(BScaledFilePath, usecols=(0, 1, 13), unpack=True, skiprows=1)
# -------- LOAD BLOS DATA --------

# -------- MANUALLY ADD ZEEMAN DATA
zeemanRa1 = 51.32
zeemanDec1 = 31.12
zeemanB1 = -27
# -------- MANUALLY ADD ZEEMAN DATA


# -------- FUNCTION DEFINITION --------
def B2RGB(b):
    """
    Takes BLOS values and assigns them a marker colour and size for use in plotting BLOS data

    :param b: The BLOS value, or list of BLOS values
    :return:  A tuple of (colour, size) corresponding to the rotation measure. Note "colour" is a tuple of (RBG,alpha)
    """
    c = []  # Marker colour
    s = []  # Marker size

    for item in b:
        if abs(item) < 1000:
            s.append(abs(item) / 2)
        if abs(item) >= 1000:
            s.append(1000 / 2)

        alpha = 1  # Optional: set the transparency
        if int(np.sign(item)) == -1:
            c.append((1, 0, 0, alpha))  # Negative rotation measures assigned red
        if int(np.sign(item)) == 1:
            c.append((0, 0, 1, alpha))  # Positive rotation measures assigned blue
        if np.sign(item) == 0:
            c.append((0, 1, 0, alpha))  # Zero-value rotation measures assigned green

    # return the list of RGBA tuples and sizes
    return c, s
# -------- FUNCTION DEFINITION. --------


# -------- READ FITS FILE --------
hdulist = fits.open(regionOfInterest.fitsFilePath)
hdu = hdulist[0]
wcs = WCS(hdu.header)
# -------- READ FITS FILE. --------

# -------- PREPARE TO PLOT BLOS POINTS --------
# ---- Convert Ra and Dec of BLOS points into pixel values of the fits file
x = []  # x pixel coordinate of BLOS
y = []  # y pixel coordinate of BLOS
for i in range(len(raDegree)):
    pixelRow, pixelColumn = wcs.wcs_world2pix(raDegree[i], decDegree[i], 0)
    x.append(pixelRow)
    y.append(pixelColumn)
# ---- Convert Ra and Dec of BLOS points into pixel values of the fits file.

color, size = B2RGB(BScaled)

n = [index+1 for index in range(len(BScaled))]
# -------- PREPARE TO PLOT BLOS POINTS. --------

# -------- PREPARE TO PLOT ZEEMAN POINTS --------
zeemanRapixel1, zeemanDecpixel1 = wcs.wcs_world2pix(zeemanRa1, zeemanDec1, 0)
# -------- PREPARE TO PLOT ZEEMAN POINTS --------

# -------- CREATE A FIGURE --------
fig = plt.figure(num=None, figsize=(12, 12), dpi=120, facecolor='w', edgecolor='k')
ax = fig.add_subplot(111, projection=wcs)

plt.title(r'$\rm{b}_{LOS}$' + ' in the Perseus Cloud\n', fontsize=12, linespacing=1)
im = plt.imshow(hdu.data, origin='lower', cmap='BrBG', interpolation='nearest')
plt.scatter(x, y, marker='o', s=size, facecolor=color, linewidth=.5)
plt.scatter(zeemanRapixel1, zeemanDecpixel1, marker="$%s$" % u'\u22A0', s=100, lw=0.5, facecolor='black')

# ---- Style the main axes and their grid
# ax.set_xlim(687, 1120)
# ax.set_ylim(275, 630)

ra = ax.coords[0]
dec = ax.coords[1]
ra.set_axislabel('\n RA (J2000)', fontsize=12, linespacing=1)
dec.set_axislabel('Dec (J2000) \n', fontsize=12, linespacing=1)

dec.set_ticks(number=10)
ra.set_ticks(number=20)
ra.display_minor_ticks(True)
dec.display_minor_ticks(True)
ra.set_minor_frequency(10)
# ra.set_ticklabel(color='w', size=10)
# dec.set_ticklabel(color='w', size=10)

ra.grid(color='black', alpha=0.5, linestyle='solid')
dec.grid(color='black', alpha=0.5, linestyle='solid')
# ---- Style the main axes and their grid.

# ---- Style the colour bar
cb = plt.colorbar(im, orientation='vertical', ticklocation='right', fraction=0.02, pad=0.145)
cb.ax.set_title(' A' + r'$_V$', linespacing=0.5, fontsize=12)
# ---- Style the colour bar.

# ---- Annotate the BLOS Points
for i, txt in enumerate(n):
    ax.annotate(txt, (x[i], y[i]), size=9, color='w')
# ---- Annotate the BLOS Points.

# ---- Style the legend
marker1 = plt.scatter([], [], s=10/2, facecolor=(1, 1, 1, 0.7))
marker2 = plt.scatter([], [], s=100/2, facecolor=(1, 1, 1, 0.7))
marker3 = plt.scatter([], [], s=500/2, facecolor=(1, 1, 1, 0.7))
marker4 = plt.scatter([], [], s=1000/2, facecolor=(1, 1, 1, 0.7))
marker5 = plt.scatter([], [], s=100, facecolor=(1, 0, 0, 0.7))
marker6 = plt.scatter([], [], s=100, facecolor=(0, 0, 1, 0.7))
marker7 = plt.scatter([], [], marker="$%s$" % u'\u22A0', s=100, lw=0.5, facecolor='black')

legend_markers = [marker1, marker2, marker4, marker5, marker6, marker7]

labels = [
    str(10)+r'$\mu G$',
    str(100)+r'$\mu G$',
    str(1000) + "+"+r'$\mu G$',
    'Away from us',
    'Towards us',
    'Zeeman' + '\n' + r'$27\mu G \, (Towards\, us)$'
    ]

legend = fig.legend(handles=legend_markers, labels=labels, scatterpoints=1, ncol=2, bbox_to_anchor=[0.433, 0.785])

frame = legend.get_frame()
frame.set_facecolor('1')
frame.set_alpha(0.4)

for label in legend.get_texts():
    label.set_fontsize('small')
# ---- Style the legend.

# ---- Display or save the figure
# plt.show()
plt.savefig(saveFilePath)
# ---- Display or save the figure.
# -------- CREATE A FIGURE. --------
