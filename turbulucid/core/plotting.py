import numpy as np
import matplotlib.pyplot as plt

__all__ = ["plot_field"]

def plot_field(field, cmap="RdBu_r", aspect="equal", **kwargs):
    """ Plot the 2d field as a contour plot."""

    pcm = plt.pcolormesh(field["X"], field["Y"], field["V"], cmap='RdBu_r', **kwargs)
    pcm.axes.set_aspect(aspect)
    plt.xlim([np.min(field["X"][:,0]), np.max(field["X"][:,-1])])
    plt.ylim([np.min(field["Y"][0,:]), np.max(field["Y"][-1,:])])
