# Import some standard python packages
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import rc
import pdb
import sys
mpl.rc('font', family='Times New Roman')
mpl.rcParams['font.size'] = 25.0

# Import coronagraph model
import coronagraph as cg

################################
# PARAMETERS
################################

# Integration time (hours)
Dt = 20.0

# Planet params
alpha = 90.     # phase angle at quadrature
Phi   = 1.      # phase function at quadrature (already included in SMART run)
Rp    = 1.0     # Earth radii
r     = 1.0     # semi-major axis (AU)

# Stellar params
Teff  = 5780.   # Sun-like Teff (K)
Rs    = 1.      # star radius in solar radii

# Planetary system params
d    = 10.     # distance to system (pc)
Nez  = 1.      # number of exo-zodis

# Plot params
plot = True
ref_lam = 0.55
saveplot = True
title = ""
ylim =  [-0.1, 0.8]
xlim =  None
tag = ""

# Save params
savefile = True

################################
# READ-IN DATA
################################

# Read-in spectrum file
fn = 'planets/earth_quadrature_radiance_refl.dat'
model = np.loadtxt(fn, skiprows=8)
lamhr = model[:,0]
radhr = model[:,1]
solhr = model[:,2]

# Calculate hi-resolution reflectivity
Ahr   = np.pi*(np.pi*radhr/solhr)

################################
# RUN CORONAGRAPH MODEL
################################

# Run coronagraph with default LUVOIR telescope (aka no keyword arguments)
lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, DtSNR = \
    cg.count_rates(Ahr, lamhr, solhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez
                   )

# Calculate background photon count rates
cb = (cz + cez + csp + cD + cR + cth)

# Convert hours to seconds
Dts = Dt * 3600.

# Calculate signal-to-noise assuming background subtraction (the "2")
SNR  = cp*Dts/np.sqrt((cp + 2*cb)*Dts)

# Calculate 1-sigma errors
sig= Cratio/SNR

# Add gaussian noise to flux ratio
spec = Cratio + np.random.randn(len(Cratio))*sig

################################
# PLOTTING
################################

if plot:

    # Create figure
    fig = plt.figure(figsize=(12,10))
    gs = gridspec.GridSpec(1,1)
    ax = plt.subplot(gs[0])

    # Set string for plot text
    if Dt > 2.0:
        timestr = "{:.0f}".format(Dt)+' hours'
    else:
        timestr = "{:.0f}".format(Dt*60)+' mins'
    plot_text = r'Distance = '+"{:.1f}".format(d)+' pc'+\
    '\n Integration time = '+timestr

    # If a reference wavelength is specified then return the SNR at that wl
    # corresponding to the integration time given
    if ref_lam:
        ireflam = (np.abs(lam - ref_lam)).argmin()
        ref_SNR = SNR[ireflam]
        plot_text = plot_text + '\n SNR = '+"{:.1f}".format(ref_SNR)+\
            ' at '+"{:.2f}".format(lam[ireflam])+r' $\mu$m'

    # Draw plot
    ax.plot(lam, Cratio*1e9, lw=2.0, color="purple", alpha=0.7, ls="steps-mid")
    ax.errorbar(lam, spec*1e9, yerr=sig*1e9, fmt='o', color='k', ms=5.0)

    # Set labels
    ax.set_ylabel(r"F$_p$/F$_s$ ($\times 10^9$)")
    ax.set_xlabel("Wavelength [$\mu$m]")
    ax.set_title(title)
    ax.text(0.99, 0.99, plot_text,\
         verticalalignment='top', horizontalalignment='right',\
         transform=ax.transAxes,\
         color='black', fontsize=20)

    # Adjust x,y limits
    if ylim is not None: ax.set_ylim(ylim)
    if xlim is not None: ax.set_xlim(xlim)

    # Save plot if requested
    if saveplot:
        plot_tag = "luvoir_demo_"+title+tag+".png"
        fig.savefig(plot_tag)
        print 'Saved: ' + plot_tag
    else:
        plt.show()

