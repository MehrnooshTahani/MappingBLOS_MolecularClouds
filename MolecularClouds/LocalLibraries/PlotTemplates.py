import matplotlib.pyplot as plt
import math

from astropy.wcs import WCS

def extinctionPlot(hdu, regionOfInterest):
    wcs = WCS(hdu.header)

    fig = plt.figure(figsize=(8, 8), dpi=120, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111, projection=wcs)

    im = plt.imshow(hdu.data, origin='lower', cmap='BrBG', interpolation='nearest')

    # ---- Style the main axes and their grid
    if not math.isnan(regionOfInterest.xmax) and not math.isnan(regionOfInterest.xmin):
        ax.set_xlim(regionOfInterest.xmin, regionOfInterest.xmax)
    if not math.isnan(regionOfInterest.ymax) and not math.isnan(regionOfInterest.ymin):
        ax.set_ylim(regionOfInterest.ymin, regionOfInterest.ymax)

    ra = ax.coords[0]
    dec = ax.coords[1]
    ra.set_major_formatter('d')
    dec.set_major_formatter('d')
    ra.set_axislabel('RA (degree)')
    dec.set_axislabel('Dec (degree)')

    dec.set_ticks(number=10)
    ra.set_ticks(number=20)
    ra.display_minor_ticks(True)
    dec.display_minor_ticks(True)
    ra.set_minor_frequency(10)

    ra.grid(color='black', alpha=0.5, linestyle='solid')
    dec.grid(color='black', alpha=0.5, linestyle='solid')
    # ---- Style the main axes and their grid.

    # ---- Style the overlay and its grid
    overlay = ax.get_coords_overlay('galactic')

    overlay[0].set_axislabel('Longitude')
    overlay[1].set_axislabel('Latitude')

    overlay[0].set_ticks(color='grey', number=20)
    overlay[1].set_ticks(color='grey', number=20)

    overlay.grid(color='grey', linestyle='solid', alpha=0.7)
    # ---- Style the overlay and its grid.

    # ---- Style the colour bar
    if regionOfInterest.fitsDataType == 'HydrogenColumnDensity':
        cb = plt.colorbar(im, ticklocation='right', fraction=0.02, pad=0.145, format='%.0e')
        cb.ax.set_title('Hydrogen Column Density', linespacing=0.5, fontsize=12)
    elif regionOfInterest.fitsDataType == 'VisualExtinction':
        cb = plt.colorbar(im, ticklocation='right', fraction=0.02, pad=0.145)
        cb.ax.set_title(' A' + r'$_V$', linespacing=0.5, fontsize=12)
    # ---- Style the colour bar.

    return fig, ax