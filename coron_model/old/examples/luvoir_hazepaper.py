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

# planet
whichplanet = 'hazyarchean'

# Integration time (hours)
Dt = 50.0
wantsnr = 20.

# Instrument params
Res = 150. # vis resolution - SLIDER
Res_UV = 20. # UV resolution - SLIDER
Res_NIR = 150. #NIR resolution - SLIDER
De_UV = 1e-4 #dark current
De_VIS = 1e-4
De_NIR = 2e-3
Re_UV = 0. #read noise
Re_VIS = 0.
Re_NIR = 0.018

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
d    = 4.5    # distance to system (pc)
dm = 4.5
Nez  = 1.      # number of exo-zodis

# Plot params
plot = True
ref_lam = 1.1
title1 = "Sun"
title2 = 'M4V star'
ylim =  [-0.05, 0.2]
xlim =  [0, 2.5]
ylim2 =  [-0.05, 15]
xlim2 =  [0, 2.5]
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
    cg.count_rates(Ahr, lamhr, solhr, alpha,  Rp, Teff, Rs, r, d, Nez, IWA=1., OWA=40., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True, diam=12.7, \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR)

lam_M, dlam_M, A_M, q_M, Cratio_M, cp_M, csp_M, cz_M, cez_M, cD_M, cR_M, cth_M, DtSNR_M = \
    cg.count_rates(Ahr, lamhr, solhrM, alpha,  Rp, Teff2, Rs, rM, dm, Nez, IWA=1., OWA=40., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True,  diam=12., \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR)

lam_7, dlam_7, A_7, q_7, Cratio_7, cp_7, csp_7, cz_7, cez_7, cD_7, cR_7, cth_7, DtSNR_7 = \
    cg.count_rates(Ahr, lamhr, solhr, alpha,  Rp, Teff, Rs, r, d, Nez, IWA=1., OWA=40., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True,  diam=7.6,  \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR)

lam_M_7, dlam_M_7, A_M_7, q_M_7, Cratio_M_7, cp_M_7, csp_M_7, cz_M_7, cez_M_7, cD_M_7, cR_M_7, cth_M_7, DtSNR_M_7 = \
    cg.count_rates(Ahr, lamhr, solhrM, alpha,  Rp, Teff2, Rs, rM, dm, Nez, IWA=1., OWA=40., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True,  diam=7.6,  \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR)

lam_5, dlam_5, A_5, q_5, Cratio_5, cp_5, csp_5, cz_5, cez_5, cD_5, cR_5, cth_5, DtSNR_5 = \
    cg.count_rates(Ahr, lamhr, solhr, alpha,  Rp, Teff, Rs, r, d, Nez, IWA=1., OWA=40., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True,  diam=5.5,  \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR)

lam_M_5, dlam_M_5, A_M_5, q_M_5, Cratio_M_5, cp_M_5, csp_M_5, cz_M_5, cez_M_5, cD_M_5, cR_M_5, cth_M_5, DtSNR_M_5 = \
    cg.count_rates(Ahr, lamhr, solhrM, alpha,  Rp, Teff2, Rs, rM, dm, Nez, IWA=1., OWA=40., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True,  diam=5.5,  \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR)

lam_4, dlam_4, A_4, q_4, Cratio_4, cp_4, csp_4, cz_4, cez_4, cD_4, cR_4, cth_4, DtSNR_4 = \
    cg.count_rates(Ahr, lamhr, solhr, alpha,  Rp, Teff, Rs, r, d, Nez, IWA=1., OWA=40., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True,  diam=4,  \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR)

lam_M_4, dlam_M_4, A_M_4, q_M_4, Cratio_M_4, cp_M_4, csp_M_4, cz_M_4, cez_M_4, cD_M_4, cR_M_4, cth_M_4, DtSNR_M_4 = \
    cg.count_rates(Ahr, lamhr, solhrM, alpha,  Rp, Teff2, Rs, rM, dm, Nez, IWA=1., OWA=40., Tsys=270.,\
                   wantsnr=wantsnr, lammin=0.2, THERMAL=True,  diam=4,  \
                   De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, \
                   Res=Res, Res_UV = Res_UV, Res_NIR = Res_NIR)

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

#baseline M dwarf 
# Calculate background photon count rates
cb_M = (cz_M + cez_M + csp_M + cD_M + cR_M + cth_M)
# Convert hours to seconds
Dts_M = Dt * 3600.
# Calculate signal-to-noise assuming background subtraction (the "2")
SNR_M  = cp_M*Dts_M/np.sqrt((cp_M + 2*cb_M)*Dts_M)
# Calculate 1-sigma errors
sig_M = Cratio_M/SNR_M
# Add gaussian noise to flux ratio
spec_M = Cratio_M + np.random.randn(len(Cratio_M))*sig_M

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
    
    # Create figure for M dwarf
    fig = plt.figure(figsize=(12,10))
    gs = gridspec.GridSpec(1,1)
    ax = plt.subplot(gs[0])

    # Set string for plot text
    if Dt > 2.0:
        timestr = "{:.0f}".format(Dt)+' hours'
    else:
        timestr = "{:.0f}".format(Dt*60)+' mins'
    plot_text = r'Distance = '+"{:.1f}".format(dm)+' pc'+\
    '\n Integration time = '+timestr

    # If a reference wavelength is specified then return the SNR at that wl
    # corresponding to the integration time given
    if ref_lam:
        ireflam = (np.abs(lam_M - ref_lam)).argmin()
        ref_SNR = SNR_M[ireflam]
        plot_text = plot_text + '\n SNR = '+"{:.1f}".format(ref_SNR)+\
            ' at '+"{:.2f}".format(lam[ireflam])+r' $\mu$m'

    # Draw plot
    ax.plot(lam_M, Cratio_M*1e9, lw=2.0, color="Salmon", alpha=0.7, ls="steps-mid")
    ax.errorbar(lam_M, spec_M*1e9, yerr=sig_M*1e9, fmt='o', color='k', ms=5.0)

    # Set labels
    ax.set_ylabel(r"F$_p$/F$_s$ ($\times 10^9$)")
    ax.set_xlabel("Wavelength [$\mu$m]")
    ax.set_title(title2)
    ax.text(0.99, 0.94, plot_text,\
         verticalalignment='top', horizontalalignment='right',\
         transform=ax.transAxes,\
            color='black', fontsize=15)

    # Adjust x,y limits
    if ylim is not None: ax.set_ylim(ylim2)
    if xlim is not None: ax.set_xlim(xlim2)

    # Save plot if requested
    if saveplot:
        plot_tag = "luvoir_demo_"+title2+tag+".png"
        fig.savefig(plot_tag)
        print 'Saved: ' + plot_tag
    #else:
    #    plt.show()
    
#~~~~~~~~~~~~~~~~~~~
    # Create Integration Time figure
    fig = plt.figure(figsize=(8.5,11))
    #gs = gridspec.GridSpec(1,1)
    #ax = plt.subplot(gs[0])
    #fig.tick_params(axis='both', which='major', labelsize=18)
    #fig.tick_params(axis='both', which='minor', labelsize=10)

   # print lam
  #  print DtSNR

    # Set string for plot text
    if Dt > 2.0:
        timestr = "{:.0f}".format(Dt)+' hours'
    else:
        timestr = "{:.0f}".format(Dt*60)+' mins'
    plot_text = r'Distance = '+"{:.1f}".format(d)+' pc'
    sun_text = 'Sun'
    gj_text = 'GJ 876'
    wantTsnr = "{:.0f}".format(wantsnr)



    # Set labels
#    fig.set_ylabel(r"Integration Time Required for SNR = "+wantTsnr, fontsize=18)
#    fig.set_xlabel("Wavelength [$\mu$m]", fontsize=18)
#    fig.set_title(title)
#    fig.text(0.99, 0.17, plot_text,\
#         verticalalignment='top', horizontalalignment='right',\
#         transform=ax.transAxes,\
#            color='black', fontsize=18)
#    fig.text(0.99, 0.12, sun_text,\
#         verticalalignment='top', horizontalalignment='right',\
#         transform=fig.transAxes,\
#            color='darkturquoise', fontsize=18)
#    fig.text(0.99, 0.07, gj_text,\
 #        verticalalignment='top', horizontalalignment='right',\
 ##        transform=fig.transAxes,\
 #           color='mediumvioletred', fontsize=18)

    sub1 = fig.add_subplot(411)
    # If a reference wavelength is specified then return the SNR at that wl
    # corresponding to the integration time given
    if ref_lam:
        ireflam = (np.abs(lam - ref_lam)).argmin()
        ref_SNR = SNR[ireflam]
        plot_text = plot_text 
    # Draw plot
    sub1.plot(lam, DtSNR, lw=3.0, color="darkturquoise", alpha=1, ls="steps-mid")
    sub1.plot(lam_M, DtSNR_M, lw=3.0, color="mediumvioletred", alpha=1, ls="steps-mid")
    sub1.axis([0.3, 2.1, 1, 10000])
    sub1.semilogy()
    sub1.set_title('15 m segmented')
    #fig.errorbar(lam, spec*1e9, yerr=sig*1e9, fmt='o', color='k', ms=5.0)
            
    sub2 = fig.add_subplot(412)
    # If a reference wavelength is specified then return the SNR at that wl
    # corresponding to the integration time given
    if ref_lam:
        ireflam = (np.abs(lam - ref_lam)).argmin()
        ref_SNR = SNR[ireflam]
        plot_text = plot_text 
    # Draw plot
    sub2.plot(lam_7, DtSNR_7, lw=3.0, color="darkturquoise", alpha=1, ls="steps-mid")
    sub2.plot(lam_M_7, DtSNR_M_7, lw=3.0, color="mediumvioletred", alpha=1, ls="steps-mid")
    sub2.axis([0.3, 2.1, 1, 10000])
    sub2.semilogy()    
    sub2.set_title('9 m segmented')
    #fig.errorbar(lam, spec*1e9, yerr=sig*1e9, fmt='o', color='k', ms=5.0)

    sub3 = fig.add_subplot(413)
    # If a reference wavelength is specified then return the SNR at that wl
    # corresponding to the integration time given
    if ref_lam:
        ireflam = (np.abs(lam - ref_lam)).argmin()
        ref_SNR = SNR[ireflam]
        plot_text = plot_text 
    # Draw plot
    sub3.plot(lam_5, DtSNR_5, lw=3.0, color="darkturquoise", alpha=1, ls="steps-mid")
    sub3.plot(lam_M_5, DtSNR_M_5, lw=3.0, color="mediumvioletred", alpha=1, ls="steps-mid")
    sub3.axis([0.3, 2.1, 1, 50000])
    sub3.semilogy()   
    sub3.set_title('6.5 m segmented')
    #fig.errorbar(lam, spec*1e9, yerr=sig*1e9, fmt='o', color='k', ms=5.0)


    sub4=fig.add_subplot(414)
    # If a reference wavelength is specified then return the SNR at that wl
    # corresponding to the integration time given
    if ref_lam:
        ireflam = (np.abs(lam - ref_lam)).argmin()
        ref_SNR = SNR[ireflam]
        plot_text = plot_text 
    # Draw plot
    sub4.plot(lam_4, DtSNR_4, lw=3.0, color="darkturquoise", alpha=1, ls="steps-mid")
    sub4.plot(lam_M_4, DtSNR_M_4, lw=3.0, color="mediumvioletred", alpha=1, ls="steps-mid")
    sub4.axis([0.3, 2.1, 1, 100000])
    sub4.semilogy()   
    sub4.set_title('4 m monolith')
    #fig.errorbar(lam, spec*1e9, yerr=sig*1e9, fmt='o', color='k', ms=5.0)
                
 
    plt.tight_layout(pad = 0.4)
            
    # Adjust x,y limits
    if ylim is not None: ax.set_ylim([0,300])
    if xlim is not None: ax.set_xlim([0.3,2.1])

    # Save plot if requested
    plot_tag = "luvoir_sun_gj876.eps"
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
