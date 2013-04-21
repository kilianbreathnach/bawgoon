import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
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
from starget import patchmk


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
        sample_list += [ abs(image[x, y] - image[x + 5, y + 5]) ]

    return np.median(np.array(sample_list))

def Initial_PSF(FWHM,double=False):

    # NB. FWHM of PSF is given in pixels.

    if not double:
        # Single Gaussian default:
        w = np.array([1.0])
        mu = np.array([[0.0,0.0]])               # centroid position in pixels
        var = (FWHM/2.35)**2.0
        cov = np.array([[[var,0.0],[0.0,var]]])  # pixels^2, covariance matrix

    else:
        # Double Gaussian alternative:
        w = np.array([0.75,0.25])
        mu = np.array([[0.0,0.0],[0.0,0.0]])
        var = (FWHM/2.35)**2.0
        cov = np.array([[[var,0.0],[0.0,var]],[[4*var,0.0],[0.0,4*var]]])

    return tractor.GaussianMixturePSF(w,mu,cov)


