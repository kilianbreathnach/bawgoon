import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import os
import logging
import sys
import numpy as np
import random as rd
from glob import glob
import pyfits
import tractor
import tractor.sdss_galaxy as sdss
from astrometry.util.util import Tan

wcsfile = "/data2/nova/BACKUP-jobs/00000046/wcs.fits"

folder = os.path.dirname(wcsfile)
pub_name = folder.split("/")[-1]    # to save on the interwebs
rdlsfiles = glob("./rdls*fits")
log = folder + "/job.log"

words = open(log).readlines()[3].split()

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

pic = picture[:,:,0].astype("float64")


# wcshead = pyfits.open(wcsfile)[0].header


#   # get the reference point of our picture, as found by Astrometry.net
#
#   refra = wcshead["CRVAL1"]
#   refdec = wcshead["CRVAL2"]
#
#   refx = wcshead["CRPIX1"]
#   refy = wcshead["CRPIX2"]
#
#   trans_mat = np.array([[wcshead["CD1_1"], wcshead["CD1_2"]],
#                         [wcshead["CD2_1"], wcshead["CD2_2"]]])
#
#   # to turn ra, dec into pixel coordinates, we'll need the inverse of the
#   # transformation matrix
#
#   imat = np.linalg.inv(trans_mat)
#
#   print imat


# get the wcs info

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

    ##   src_radec = np.asarray(rd.choice(cat_srcs))
    #
    #    '''get difference in ra, dec from reference point'''
    #    radec_diff = src_radec - np.array([refra, refdec])
    #
    #    '''now turn that difference into a pixel difference'''
    #    pix_diff = np.dot(imat, radec_diff)
    #
    #    '''now add the pixel difference to the pixel reference to get X,Y coords'''
    #    src_xy = pix_diff + np.array([refx, refy])
    #
    #    centre = np.rint(src_xy).astype("int")
    #
    #    starr = pic[centre[0] - 3:centre[0] + 4, centre[1] - 3:centre[1] + 4]

    ra, dec = src

    print "RA: ", ra
    print "DEC: ", dec

    src_xy = np.asarray(wcs.wcs.radec2pixelxy(ra, dec))
    print src_xy

    centre = np.rint(src_xy).astype("int")
    print "centre:"
    print centre[0]
    print centre[1]
    starr = pic[centre[0] - 12:centre[0] + 13, centre[1] -
            12:centre[1] + 13]

    print starr

    # and let's have a go at plotting so

    fig = plt.figure()

    ax = fig.add_subplot(111)
    try:
        ax.imshow(starr, cmap=plt.cm.gray, interpolation="nearest")
    except:
        continue

    fig.savefig("/home/kilian/public_html/exp_req/test/7x7_{0}.png".format(run))
