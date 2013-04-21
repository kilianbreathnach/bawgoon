import starget as sg


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




f = open("stars.dat", mode="a")

for img in range(1):

    #    img =

    try:
        fits = pyfits.open(goodfile)
        picture = fits[0].data
        if type(picture) != np.ndarray:
            picture = fits[1].data
    except (IOError):
        picture = mpimg.imread(goodfile)[::-1]

    pic = picture[:,:,0].astype("float64")


    for star in range(1):


