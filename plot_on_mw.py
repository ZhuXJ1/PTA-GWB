
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg
from astropy import units as u
from astropy.coordinates import SkyCoord

# inputs:
d_sun = 8.5  # kpc
mw_png = 'your_image'
fs = 12  # fontsize in your plot

### plot nord 42 emi on the artist milky way image ###
fig, ax = plt.subplots()
# read png img: http://www.scipy-lectures.org/advanced/image_processing/
img = mpimg.imread(mw_png)
plt.imshow(img)

# the location of the Sun and the Galactic center in the unit of pixel
x_sun = 2800
y_gc = 2800
x_gc = 2800
y_sun = 3870
r_mw = 20  # radius of the Milky Way disk

pix_per_kpc = (y_sun - y_gc) / d_sun  # for conversion from kpc to pixel

def kpc_to_pix(x_kpc, y_kpc):
    """
    convert position in unit of kpc to pixel
    The center of the coordiate on the image is in the location of the Sun
    :param x_kpc: the x value of your pulsars in kpc, a float or a numpy array
    :param y_kpc: the y value of your pulsars in kpc, a float or a numpy array
    :return: 
    """
    x_pix = x_sun + x_kpc * pix_per_kpc
    y_pix = y_sun - y_kpc * pix_per_kpc
    return x_pix, y_pix


# calculate positions of the starting and ending points
x_hii_pix, y_hii_pix = kpc_to_pix(x_hii_kpc, y_hii_kpc)
x_los_pix, y_los_pix = kpc_to_pix(x_los_kpc, y_los_kpc)

plt.plot(x_sun, y_sun, 'w+', markersize=8)  # plot the Sun
plt.plot(x_gc, y_gc, 'k+', markersize=8)  # plot the Galactic center
ax.text(x_gc - 200, y_gc - 150, 'GC', fontsize=fs)
ax.text(x_sun + 200, y_sun + 100, 'Sun', fontsize=fs, color='w')

# this is the lines I plotted
data_back_nord = [(x_hii_pix, x_los_pix), (y_hii_pix, y_los_pix)]


# below are for plotting colorbars
my_cmap = 'rainbow'  # backup color maps: gist_heat, rainbow, cubehelix
a = plt.get_cmap(my_cmap)  # get color values
emax = np.max(emib_nord)
emin = np.min(emib_nord)
ax.set_color_cycle([a((i - emin) / (emax - emin)) for i in emib_nord])  # make i range in (0,1)
plt.plot(*data_back_nord, linewidth=0.6)
# plot color bar
sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.Normalize(vmin=emin, vmax=emax))
sm._A = []
cbar = plt.colorbar(sm, format='%.1f')
cbar.ax.set_ylabel('Emissivity (K pc$^{-1}$)', fontsize=fs)
cbar.ax.tick_params(labelsize=fs)

# plot the white dashed circle, this is the edge of the Milky Way disk
g_edge = plt.Circle((x_gc, y_gc), r_mw * pix_per_kpc, color='white', fill=False, linestyle='dashed')
fig.gca().add_artist(g_edge)
fig = plt.gcf()
fig.gca().add_artist(g_edge)
plt.axis('off')
plt.savefig(path_measure + 'fig/emi_on_mw_img_nord.png', dpi=300, bbox_inches="tight")




def RA_DEC_to_GLong_GLat(RA, DEC):
     """change the coordinates from (RA, DEC) in degree to (GLong, GLat) in degree"""
     # python 2.7
     """
     a = ICRS(ra=RA, dec=DEC, unit=(u.degree, u.degree))
     b = a.transform_to(Galactic)
     """
     # python 3.6
     c_icrs = SkyCoord(ra=RA * u.degree, dec=DEC * u.degree, frame='icrs')
     c_galactic = c_icrs.transform_to('galactic')  # c_icrs.galactic does the same thing
     return float(c_galactic.l.degree), float(c_galactic.b.degree)