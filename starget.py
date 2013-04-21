#########################################################################
# This module processes the desired star. It contains functions to find #
# the star given a data image and a catalogue reference. The main       #
# function is patchmk, which generates a 7x7 (?) patch image centred on #
# the star's location in the data image as proposed by the catalogue    #
# This patch will then be tested for compatibility of image and         #
# catalogue.                                                            #
#########################################################################

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

def patschmk(img, wcs, x, y):

    return patschi
