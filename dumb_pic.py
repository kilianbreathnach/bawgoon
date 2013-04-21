import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import os
import logging
import sys
import numpy as np
import random as rd
import collections as colxns
from glob import glob
import pyfits
import tractor
import tractor.sdss_galaxy as sdss
from astrometry.util.util import Tan
from scipy.optimize import curve_fit as cf
from scipy.optimize import leastsq as lsq
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


MAGIC_NO = 1.0484
'''For converting MAD to standard deviation.'''


def MAD(image):
    '''
    Get Median Absolute Deviation of image.
    image should be an np array.
    '''
    row, col = image.shape
    sample_list = []
    for i in range(50):
        x, y = rd.randrange(0, row - 6), rd.randrange(0, col - 5)
        sample_list += [abs(image[x, y] - image[x + 5, y + 5])]

    return np.median(np.array(sample_list))


def residuals(params, func, x, data):
    res = data - func(params, x)
    return res.reshape(np.prod(res.shape))


def sky_mod(params, base):
    """
    - plane that defines the sky background with tilt
    - v is vertical tilt and h is horizontal tilt
    - angles are positive for brighter top/left respectively
    - x is a 625 element long array of ones
    """
    z0, mh, mv = params
    sky = z0 * base.reshape(25, 25)

    for i in range(25):
        for j in range(25):
            #dy = np.abs(j - 12)
            sky[j][i] += i * mh + j * mv  #dx * mh + dy * mv

    return sky


def psf_mod(params, base):
    '''
    Model a background sky with tilt and a bivariate gaussian PSF
    - base should be a 625-element array of ones
    - means of gaussians both taken as zero (centred on central pixel)
    - see sky_mod for explanation of sky parameters
    - w is weight which gaussian has compared to sky
    - varx, vary are x, y gaussian variances
    - cor is correlation of variances
    '''
    z0, mh, mv, w, varx, vary, cor = params
    sky = base.reshape(25, 25) * z0
    q = 1. - cor ** 2
    gauss = np.ones((25, 25)) * (w / (2 * np.pi * np.sqrt(varx * vary * q)))

    for i in range(25):
        ipart = np.exp((-0.5 / q) * ((i - 12) ** 2 / varx))
        for j in range(25):
            sky[j][i] += i * mh + j * mv
            gauss[j][i] *= ipart * np.exp((-0.5 / q) * (((j - 12)**2 / vary) -
                                                ((2.*cor*(i - 12)*(j - 12)) /
                                                 np.sqrt(varx * vary))))

    return sky + gauss


def roam_mod(params, base):
    '''
    Model a background sky with tilt and a bivariate gaussian PSF
    - base should be a 625-element array of ones
    - means of gaussians both taken as zero (centred on central pixel)
    - see sky_mod for explanation of sky parameters
    - w is weight which gaussian has compared to sky
    - varx, vary are x, y gaussian variances
    - cor is correlation of variances
    '''
    z0, mh, mv, w, varx, vary, cor, x, y = params
    sky = base.reshape(25, 25) * z0
    q = 1. - cor ** 2
    gauss = np.ones((25, 25)) * (w / (2 * np.pi * np.sqrt(varx * vary * q)))

    for i in range(25):
        ipart = np.exp((-0.5 / q) * ((i - x)**2 / varx))
        for j in range(25):
            sky[j][i] += i * mh + j * mv
            gauss[j][i] *= ipart * np.exp((-0.5 / q) * (((j - y)**2 / vary) -
                                                ((2.*cor*(i - x)*(j - y)) /
                                                 np.sqrt(varx * vary))))

    return sky + gauss


def chi_img(data, model, sig):

    return (data - model) / sig



wcsfile = "/data2/nova/BACKUP-jobs/00000046/wcs.fits"

folder = os.path.dirname(wcsfile)
pub_name = folder.split("/")[-1]    # to save on the interwebs
rdlsfiles = glob("./rdls*fits")
log = folder + "/job.log"

words = open(log).readlines()[3].split()
print words

imagef_ind = words.index("--image") + 1    # wrong filename for the image
badfile = words[imagef_ind]
getpath = badfile.split("/")[-4:]
goodfile = "/home/nova/BACKUP-data/{0}/{1}/{2}/{3}".format(*tuple(getpath))
'''Now we have the right filename for the original.'''

try:
    fits = pyfits.open(goodfile)
    picture = fits[0].data
    if type(picture) != np.ndarray:
        picture = fits[1].data
except (IOError):
    picture = mpimg.imread(goodfile)[::-1]


# let's just work on one of the colours:

pic = picture[:, :, 0].astype("float64")

stdev = MAGIC_NO * MAD(pic)

print "dimensions are ", pic.shape

ext = 0
wcs = tractor.FitsWcs(Tan(wcsfile, ext))


# pick the ra,dec of a star to check

cat_srcs = []

for rdfile in rdlsfiles:
    rdls = pyfits.open(rdfile)
    radec = rdls[1].data
    for ra in radec:
        cat_srcs.append(ra)

for run, src in enumerate(cat_srcs):
    print "###################################"
    print "######## RUN NO. ", run, "#########"
    print "###################################"

    ra, dec = src
    t_radec = tractor.RaDecPos(ra, dec)
    print "RaDecPos obj: ", t_radec

    print wcs.positionToPixel(t_radec)

    src_xy = np.asarray(wcs.wcs.radec2pixelxy(ra, dec))
    print src_xy

    centre = np.rint(src_xy).astype("int")
    print "centre:"
    print centre[0]
    print centre[1]

     # let's just switch the centre values and see...
    starr = pic[centre[1] - 12:centre[1] + 13, centre[0] - 12:centre[0] + 13]

    print "starr is  ", starr.shape


    if starr.shape == (25, 25):

        # now try to fit to sky model

        med = np.median(starr)
        peak = np.amax(starr)
        print "median of starr: ", med
        print "max of starr: ", peak
        """guess where gaussian is centred"""
        x0, y0 = (float(pos[0]) for pos in np.where(starr == peak))

        xdata = np.ones(625)

        """
        Define plotstuff object
        - 0th term is sky only model
        - 1st term is sky and psf
        - 2nd is above plus roaming psf
        """

        skyp = colxns.namedtuple('skyp', 'z0 mh mv')
        skypsf = colxns.namedtuple("skypsf", "z0 mh mv w varx vary cor")
        roamp = colxns.namedtuple("roamp", "z0 mh mv w varx vary cor x y")

        kwargs1 = dict(z0=med, mh=0., mv=0.)
        kwargs2 = dict(kwargs1, w=peak, varx=4., vary=4., cor=0.)
        kwargs3 = dict(kwargs2, x=x0, y=y0)

        class RunParams():
            pass

        stuff = [RunParams(), RunParams(), RunParams()]

        stuff[0].tup, stuff[0].kws, stuff[0].fun = skyp, kwargs1, sky_mod
        stuff[1].tup, stuff[1].kws, stuff[1].fun = skypsf, kwargs2, psf_mod
        stuff[2].tup, stuff[2].kws, stuff[2].fun = roamp, kwargs3, roam_mod


        # and start up the figure for plotting
        fig = plt.figure()
        ax = []

        for no, plst in enumerate(stuff):
            """
            and here we do the fitting and plotting, one row for each model
            """
            # first the real image
            plt.subplot(3, 3, no * 3 + 1)
            plt.imshow(starr, cmap=plt.cm.gray, interpolation="nearest")

            # now fit and plot the model
            param_guess = plst.tup(**plst.kws)

            result = lsq(residuals, param_guess, args=(plst.fun, xdata, starr))

            solve_params = result[0]
            soln = plst.tup(*solve_params)
            print "optimal params are: ", soln
            bestsky = plst.fun(soln, xdata)

            plt.subplot(3, 3, no * 3 + 2)
            plt.imshow(bestsky, cmap=plt.cm.gray, interpolation="nearest")

            # and finally the chi image
            chimg = chi_img(starr, bestsky, stdev)
            plt.subplot(3, 3, no * 3 + 3)
            plt.imshow(chimg, cmap=plt.cm.gray, interpolation="nearest")


        plt.tight_layout()

        fig.savefig("/home/kilian/public_html/exp_req/switch-test/{0}.png".format(run))

    else:
        print "out of range"
        continue
