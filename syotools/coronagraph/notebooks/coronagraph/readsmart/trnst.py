import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
mpl.rc('font',family='Times New Roman')
mpl.rcParams['font.size'] = 25.0

def trnst(path, getdata=False, plot=True, save=False, ptype="depth", xran=None, xlog=False, ptitle=""):
    """Reads-in SMART transit output files. Plots, saves, and returns arrays if requested.

    Parameters
    ----------
    path : str
        Location of .trnst file
    getdata : bool (optional)
        Flag for whether or not to return data
        See Returns below
    plot : bool (optional)
        Flag for whether or not to plot data
    save : bool (optional)
        Flag for whether or not to save plot
    pytpe : str (optional)
        Type of plot: 'depth' or 'absrad'
    xran : list (optional)
        x-axis range for plot
    xlog : bool (optional)
        Flag for logarithmic x-axis in plot
    ptitle : str (optional)
        String for plot title

    Returns
    -------
    wl : numpy.ndarray
        Wavelength grid for transmission spectrum (um)
    wn : numpy.ndarray
        Wavenumber grid for transmission spectrum (cm^-1)
    zeff : numpy.ndarray
        Effective absorbing radius of planet (km)
    tdepth : numpy.ndarray
        Transit depth, ie change in flux, ie (Rp/Rs)^2

    """


    # Read in .trnst file
    data = np.genfromtxt(path)
    # Split into arrays
    wl = data[:,0]
    wn = data[:,1]
    zeff = data[:,2]
    tdepth = data[:,3]

    # Plot if requested
    if plot:
        fig = plt.figure(figsize=(12,10))
        gs = gridspec.GridSpec(1,1)
        ax0 = plt.subplot(gs[0])
        ax0.set_title(ptitle)
        # Plot requested type
        if ptype == "depth":
            ax0.plot(wl, tdepth, c='k', lw=2.0)
            ax0.set_ylabel(r"$(R_p / R_{\star})^2$")
        elif ptype == "absrad":
            ax0.plot(wl, zeff, c='k', lw=2.0)
            ax0.set_ylabel(r"Absorbing Radius [km]")
        else:
            print("Error: Invalid ptype. Use 'depth' or 'absrad'.")
        ax0.set_xlabel(r"Wavelength [$\mu$m]")
        if xlog: ax0.semilogx()
        if xran is not None: ax0.set_xlim(xran)
        # Save plot if requested
        if save:
            savetag = 'plot_trnst_'+ptype+'_'+ptitle+'.pdf'
            fig.savefig(savetag)
            print 'Saved plot: '+savetag
        fig.show();

    # Return data if requested
    if getdata:
        return wl, wn, zeff, tdepth
    else:
        return
