import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import numpy as np

dippy = 60.

def compare(data, model, res, chi, filename):

    head = np.percentile(data, 99.9)
    toe = np.percentile(data, 0.1)
    chirange = 5. # head*0.05 # np.percentile(data, 9.3)
    resrange = (head-toe)*0.5 # np.percentile(data, 50.0)

    width, height = np.shape(data)
    figw, figh = (width/dippy)*2., (height/dippy)*2.

    fig = plt.figure(figsize=(figw, figh))
#    fig = plt.figure(dpi=dippy, figsize=(figw, figh))

    gs = gridspec.GridSpec(2, 2)

#    axr = fig.add_subplot(221)
    axr = plt.subplot(gs[0])
    axr.set_axis_off()
    axr.set_title(r"Data", fontsize=40)

    caxr = axr.imshow(data, cmap=plt.cm.gray, vmin=toe, vmax=head)
    divr = make_axes_locatable(plt.gca())
    cbarr = divr.append_axes("right", "5%", pad="3%")
    fig.colorbar(caxr, cax=cbarr)

#    axf = fig.add_subplot(222)
    axf = plt.subplot(gs[1])
    axf.set_axis_off()
    axf.set_title(r"Model Image", fontsize=40)

    caxf = axf.imshow(model, cmap=plt.cm.gray, vmin=toe, vmax=head)
    divf = make_axes_locatable(plt.gca())
    cbarf = divf.append_axes("right", "5%", pad="3%")
    fig.colorbar(caxf, cax=cbarf)

#    axres = fig.add_subplot(223)
    axres = plt.subplot(gs[2])
    axres.set_axis_off()
    axres.set_title(r"Residual Image",fontsize=40)

    caxres = axres.imshow(res, cmap=plt.cm.gray, vmin=-resrange, vmax=resrange)
    divres = make_axes_locatable(plt.gca())
    cbarres = divres.append_axes("right", "5%", pad="3%")
    fig.colorbar(caxres, cax=cbarres)

#    axchi = fig.add_subplot(224)
    axchi = plt.subplot(gs[3])
    axchi.set_axis_off()
    axchi.set_title(r"Chi Image", fontsize=40)

    caxchi = axchi.imshow(chi, cmap=plt.cm.gray, vmin=-chirange, vmax=chirange)
    divchi = make_axes_locatable(plt.gca())
    cbarchi = divchi.append_axes("right", "5%", pad="3%")
    fig.colorbar(caxchi, cax=cbarchi)

#    fig.subplots_adjust(hspace=0.001)
    fig.tight_layout()

#    fig.savefig(filename, dpi=60., bbox_inches='tight')
    fig.savefig(filename)

