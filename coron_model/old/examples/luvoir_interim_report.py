# Import some standard python packages
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import rc
import pdb
import sys
mpl.rc('font', family='Helvetica')
mpl.rcParams['font.size'] = 25.0

# Import coronagraph model
import coronagraph as cg

################################
# PARAMETERS
################################

# planet
whichplanet = 'earth'

# Integration time (hours)
Dt = 50.0
wantsnr = 10.

# Instrument params
Res = 150. # vis resolution - SLIDER
Res_UV = 10. # UV resolution - SLIDER
Res_NIR = 150. #NIR resolution - SLIDER
De_UV = 1e-4 #dark current
De_VIS = 1e-4
De_NIR = 2e-3
Re_UV = 0. #read noise
Re_VIS = 0.
Re_NIR = 0.018
Tput = 0.05 #throughput 

# Planet params
alpha = 90.     # phase angle at quadrature
Phi   = 1.      # phase function at quadrature (already included in SMART run)
Rp    = 1.0     # Earth radii
r     = 1.0     # semi-major axis (AU)
rM   = 0.12    # semi-major axis of M dwarf planet (AU)

# Stellar params
Teff  = 5780.   # Sun-like Teff (K)
Teff2 = 3130.
Rs    = 1.      # star radius in solar radii

# Planetary system params
d    = 5.    # distance to system (pc)
dm = 4.5
Nez  = 1.      # number of exo-zodis

# Plot params
plot = True
ref_lam = 1.1
title1 = "Sun"
title2 = 'M4V star'
ylim =  [-0.05, 0.2]
xlim =  [0.41, 1.7]
ylim2 =  [-0.05, 15]
xlim2 =  [0.41, 1.7]
tag = ""

# Save params
savefile = True
saveplot = True


################################
# READ-IN DATA
################################

# Read-in spectrum file
if whichplanet == 'earth':
    fn = '../coronagraph/planets/earth_quadrature_radiance_refl.dat'
    model = np.loadtxt(fn, skiprows=8)
    lamhr = model[:,0]
    radhr = model[:,1]
    solhr = model[:,2]

    # Calculate hi-resolution reflectivity
    Ahr   = np.pi*(np.pi*radhr/solhr)

if whichplanet == 'hazyarchean':
    fn = '../coronagraph/planets/Hazy_ArcheanEarth_geo_albedo.txt'
    model = np.loadtxt(fn, skiprows=8)
    lamhr = model[:,0]
    Ahr = model[:,1]
    solhr =  cg.noise_routines.Fstar(lamhr, Teff, Rs, r, AU=True) #sun
    solhrM =  cg.noise_routines.Fstar(lamhr, 3130., 0.3761, 0.12, AU=True) #gj876
    



    
################################
# RUN CORONAGRAPH MODEL
################################

# Run coronagraph with default LUVOIR telescope (aka no keyword arguments)
#lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, DtSNR = \
#    cg.count_rates(Ahr, lamhr, solhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez\
#                   )

lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, DtSNR = \
    cg.count_rates(Ahr, lamhr, solhr, alpha,  Rp, Teff, Rs, r, d, Nez, IWA=2., OWA=64., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True, diam=15., \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR, Tput=Tput)

lam_9, dlam_9, A_9, q_9, Cratio_9, cp_9, csp_9, cz_9, cez_9, cD_9, cR_9, cth_9, DtSNR_9 = \
    cg.count_rates(Ahr, lamhr, solhr, alpha,  Rp, Teff, Rs, r, d, Nez, IWA=2., OWA=64., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True,  diam=9.,  \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR, Tput=Tput)


lam_65, dlam_65, A_65, q_65, Cratio_65, cp_65, csp_65, cz_65, cez_65, cD_65, cR_65, cth_65, DtSNR_65 = \
    cg.count_rates(Ahr, lamhr, solhr, alpha,  Rp, Teff, Rs, r, d, Nez, IWA=2., OWA=64., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True,  diam=6.5,  \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR, Tput=Tput)

lam_4, dlam_4, A_4, q_4, Cratio_4, cp_4, csp_4, cz_4, cez_4, cD_4, cR_4, cth_4, DtSNR_4 = \
    cg.count_rates(Ahr, lamhr, solhr, alpha,  Rp, Teff, Rs, r, d, Nez, IWA=2., OWA=64., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True,  diam=4,  \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR, Tput=Tput)


#baseline solar
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
    ax.plot(lam, Cratio*1e9, lw=2.0, color="green", alpha=0.7, ls="steps-mid")
    ax.errorbar(lam, spec*1e9, yerr=sig*1e9, fmt='o', color='k', ms=5.0)

    # Set labels
    ax.set_ylabel(r"F$_p$/F$_s$ ($\times 10^9$)")
    ax.set_xlabel("Wavelength [$\mu$m]")
    ax.set_title(title1)
    ax.text(0.01, 0.94, plot_text,\
         verticalalignment='top', horizontalalignment='left',\
         transform=ax.transAxes,\
            color='black', fontsize=15)

    # Adjust x,y limits
    if ylim is not None: ax.set_ylim(ylim)
    if xlim is not None: ax.set_xlim(xlim)

    # Save plot if requested
    if saveplot:
        plot_tag = "luvoir_demo_"+title1+tag+".png"
        fig.savefig(plot_tag)
        print 'Saved: ' + plot_tag
    #else:
    #    plt.show()

    #--------
    
    
#~~~~~~~~~~~~~~~~~~~
   # Create Integration Time figure
    fig = plt.figure(figsize=(6,3.6))
    #gs = gridspec.GridSpec(1,1)
    ax = plt.subplot(gs[0])
    ax.axis('off')
    #fig.tick_params(axis='both', which='major', labelsize=18)
    #fig.tick_params(axis='both', which='minor', labelsize=10)

   # print lam
  #  print DtSNR

    # Set string for plot text
    if Dt > 2.0:
        timestr = ''
    else:
        timestr = ''
    _15m_text= '15 m'
    _9m_text = '9 m'
    _65m_text = '6.5 m'
    _5m_text = '4 m'



    # Set labels
  #  ax.set_ylabel("Integration Time Required for SNR = "+wantTsnr +" hours" , fontsize=6)
#    fig.set_xlabel("Wavelength [$\mu$m]", fontsize=18)
#    fig.set_title(title)
    fig.text(0.92, 0.45, _15m_text,\
         verticalalignment='top', horizontalalignment='right',\
         transform=ax.transAxes,\
             color='darkturquoise', fontsize=10)
    fig.text(0.92, 0.4, _9m_text,\
         verticalalignment='top', horizontalalignment='right',\
         transform=ax.transAxes,\
             color='blue', fontsize=10)    
    fig.text(0.92, 0.35, _65m_text,\
         verticalalignment='top', horizontalalignment='right',\
         transform=ax.transAxes,\
             color='green', fontsize=10)
    fig.text(0.92, 0.3, _5m_text,\
         verticalalignment='top', horizontalalignment='right',\
         transform=ax.transAxes,\
             color='mediumvioletred', fontsize=10)

    sub1 = fig.add_subplot(111)
    # If a reference wavelength is specified then return the SNR at that wl
    # corresponding to the integration time given
    if ref_lam:
        ireflam = (np.abs(lam - ref_lam)).argmin()
        ref_SNR = SNR[ireflam]
        plot_text = plot_text 
    # Draw plot
    sub1.plot(lam, DtSNR, lw=1.5, color="darkturquoise", alpha=1, ls="steps-mid")
    sub1.plot(lam_9, DtSNR_9, lw=1.5, color="blue", alpha=1, ls="steps-mid")
    sub1.plot(lam_65, DtSNR_65, lw=1.5, color="green", alpha=1, ls="steps-mid")
    sub1.plot(lam_4, DtSNR_4, lw=1.5, color="mediumvioletred", alpha=1, ls="steps-mid")

   # sub1.plot(lam_M, DtSNR_M, lw=2.0, color="mediumvioletred", alpha=1, ls="steps-mid")
    sub1.axis([0.4, 1.8, 1, 100000])
    sub1.semilogy()
    sub1.tick_params(labelsize=8)
    sub1.set_title('', fontsize=8)
   # sub1.text(0.3, 100, "Integration Time Required for SNR = 10 [hours]", rotation='vertical',  horizontalalignment='center', verticalalignment='center', fontsize=12)
  #  sub1.axvline(1.64, color='red', linestyle='--')
  #  sub1.axvline(0.82, color='green', linestyle='--')
  #  sub1.axvline(0.54, color='blue', linestyle='--')
    sub1.set_xlabel("Wavelength [$\mu$m]", fontsize=8)
    sub1.set_ylabel("Integration Time Required for SNR = 10 [hours]", fontsize=8)

  
    plt.tight_layout(pad = 0.)
            
    # Adjust x,y limits
#    if ylim is not None: ax.set_ylim([1e0,1e5])
#    if xlim is not None: ax.set_xlim([0.3,2.1])

    # Save plot if requested
    plot_tag = "luvoir_integration_interim_report.eps"
    fig.savefig(plot_tag, format='eps', dpi=1000)
    print 'Saved: ' + plot_tag

#    plt.show()
################################
# SAVING
################################

# Save Synthetic data file (wavelength, albedo, error) if requested
if savefile:
    data_tag = 'luvoir_demo_'+tag+'.txt'
    y_sav = np.array([lam,spec,sig])
    np.savetxt(data_tag, y_sav.T)
    print 'Saved: ' + data_tag

sys.exit()
