
import numpy as np

def atm(path, skiprows=1, getdata=False, plot=True, ptitle="", save=False):

    '''Reads-in *.atm file. Plots and returns arrays upon request.
     
    Parameters
    ----------
    file : string
        Location/name of *.atm file.
        Ex: 'icrccm_62.atm'
    skiprows : int
        Number of rows at top of atm file before data begin.
    getdata : bool
        If True (default: False) then function returns arrays listed below.
    plot : bool
        If True (default: True) then function plots atm data.
    ptitle : string
        Optional name for generated plot

    Returns
    ----------
    P : 1-D array of floats
        Pressure grid [Pa]
    T : 1-D array of floats
        Temperature grid [K]
    abund_profiles : 2-D array of floats
        Mixing ratio grid 
        ex: >>> N_absorbers = len(abund_profile[0,:]) 
    gasform : 1-D array of strings
        List of gas formulas 

    Revision History
    ----------
    Written by J. Lustig-Yaeger Sept., 2015

    '''   	

    # Note Pressure units are assumed to be in Pascals
    
    from gas_info import gas_info
    gases = gas_info()
    
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib import gridspec
    from matplotlib import rc
    mpl.rcParams['font.size'] = 15.0

    if not getdata and not plot:
        print "Error in read_atm: Set getdata=True or plot=True, otherwise this function does nothing."
        return
    
    # Get atm data
    atm_data = np.genfromtxt(path, skip_header=skiprows-1)
    gas_codes = atm_data[0,2:]
    abund_profiles = atm_data[1:,2:]
    P = atm_data[1:,0]
    T = atm_data[1:,1]
    
    # Convet gas codes to formulae
    gasform = []
    for i in gas_codes:
        gasform.append(gases['Formula'][int(i)])
        
    # Plot...
    #'''
    if plot:
        fig = plt.figure(figsize=(18,9))
        gs = gridspec.GridSpec(1,2)
        ax0 = plt.subplot(gs[0])
        for i in range(abund_profiles[0].size):
            ax0.plot(abund_profiles[:,i], P, lw=2.0, label=tex_molecule(gasform[i]))
        ax0.set_xlabel(r'Mixing Ratio')
        ax0.set_ylabel(r'P [Pa]')
        ax0.invert_yaxis()
        ax1 = plt.subplot(gs[1])
        ax1.invert_yaxis()
        ax1.yaxis.tick_right() 
        ax1.yaxis.set_label_position("right")
        #ax1.axes.get_yaxis().set_visible(False)
        tp = ax1.plot(T,P*1e-5,lw=2.0, label='T', color='k')
        ax1.set_xlabel(r'T [K]')
        ax1.set_ylabel(r'P [bars]')
        ax0.set_xlim([np.min(abund_profiles)/10., np.max(abund_profiles)*10.])
        ax0.semilogx()
        ax0.semilogy()
        ax1.semilogy()
        leg = ax0.legend(loc=1)
        leg.get_frame().set_alpha(0.8)
        plt.tight_layout()
        plt.suptitle(ptitle)
	if save:
	    savetag = 'plot_atm_'+ptitle+'.pdf'
            fig.savefig(savetag)
	    print 'Saved plot: '+savetag
        plt.show();
    #'''
    
    
    if getdata:
        return P, T, abund_profiles, gasform
    else:
        return

#################################################################################

def rad(path, plot=False, getdata=True, ptype='TOA', ptitle='', xran=[1,5], save=False):

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib import gridspec
    from matplotlib import rc
    mpl.rcParams['font.size'] = 15.0
    
    rad_data = np.genfromtxt(path, skip_header=0)
    wl = rad_data[:,0]
    wno = rad_data[:,1]
    solar_spec = rad_data[:,2]
    TOA_flux = rad_data[:,3] #W/m**2/um
    rad_streams = rad_data[:,4::]
    
    plot_rad=False
    
    if ptype == 'TOA':
        y = TOA_flux
        ytitle = r'TOA Flux [W m$^{-2} \mu $m$^{-1}$]'
    elif ptype == 'reflect':
        y = TOA_flux / solar_spec
        ytitle = r'Reflectance'
    elif ptype == 'solar':
        y = solar_spec
        ytitle = r'Solar Flux [W m$^{-2} \mu $m$^{-1}$]' 
    elif ptype == 'rad':
        y = rad_streams
        plot = False
        plot_rad = True
    else:
        plot = False
        print "Invalid ptype. Must choose betweem 'TOA', 'reflect', 'solar', and 'rad'."
    
    if plot:
        x = wl[(wl > xran[0]) & (wl < xran[1])]
        y = y[(x > xran[0]) & (x < xran[1])]
        fig = plt.figure(figsize=(10,7))
        gs = gridspec.GridSpec(1,1)
        ax0 = plt.subplot(gs[0])
        ax0.plot(x, y, lw=2.0, label=ptitle)
        ax0.set_xlabel(r'$\lambda$ $[\mu m]$')
        ax0.set_ylabel(ytitle)
        ax0.set_title('')
        #ax0.semilogx()
        ax0.set_xlim([np.min(x), np.max(x)])
        ax0.set_ylim([np.min(y), np.max(y)])
	if save:
	    savetag = 'plot_rad_'+ptype+'_'+ptitle+'.pdf'
            fig.savefig(savetag)
            print 'Saved plot: '+savetag
	fig.show();
        
    if plot_rad:
        x = wl[(wl > xran[0]) & (wl < xran[1])]
        y = y[(x > xran[0]) & (x < xran[1])]
        num = len(rad_streams[0,:])
        fig = plt.figure(figsize=(10,7))
        gs = gridspec.GridSpec(1,1)
        ax0 = plt.subplot(gs[0])
        for i in range(num):
            ax0.plot(x, y[:,i], lw=2.0, label='Stream '+str(i+1))
        ax0.set_xlabel(r'$\lambda$ $[\mu m]$')
        ax0.set_ylabel(r'Radiance [W m$^{-2} \mu $m$^{-1}$ sr$^{-1}$]')
        leg = ax0.legend(loc=1)
        leg.get_frame().set_alpha(0.8)
	if save:
            savetag = 'plot_rad_'+ptype+'_'+ptitle+'.pdf'
            fig.savefig(savetag)
            print 'Saved plot: '+savetag
	fig.show();
        
    if getdata:
        return wl, wno, solar_spec, TOA_flux, rad_streams
    else:
        return

#################################################################################

def tran(path):
    data = np.genfromtxt(path, skip_header=3)
    wl = data[:,0]
    fluxratio = data[:,1]
    absorbing_radius = data[:,2]
    return wl, fluxratio, absorbing_radius

#################################################################################

def hrt(path, Nstreams=8, plot=False, save=False, ptitle=''):
    
    # Read first line separately since it's missing a value
    f = open(path, 'r')
    for i in range(11):
        line = f.readline()
    l = np.array(map(float, line.split()))
    f.close()

    # Read the rest, adding first value
    data = np.genfromtxt(path, skip_header=11)
    pressure = np.hstack([l[0],data[:,0]]) 
    temperature = np.hstack([l[1],data[:,1]]) 
    altitude = np.hstack([l[2],data[:,2]])  
    solarQ = data[:,3] 
    thermalQ = data[:,4] 
    dir_sol_flx = np.hstack([l[3],data[:,5]]) 
    dn_sol_flx = np.hstack([l[4],data[:,6]]) 
    up_sol_flx = np.hstack([l[5],data[:,7]]) 
    dn_th_flx = np.hstack([l[6],data[:,8]]) 
    up_th_flx = np.hstack([l[7],data[:,9]])  
    p_flux = np.hstack([l[8],data[:,10]])     
    t_flux = np.hstack([l[9],data[:,11]]) 
    
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        from matplotlib import gridspec
        from matplotlib import rc
        mpl.rcParams['font.size'] = 15.0
        fig = plt.figure(figsize=(12,8))
        gs = gridspec.GridSpec(1,1) 
        ax0 = plt.subplot(gs[0])
        ax0.plot(dir_sol_flx,altitude,c='orange', lw=2.0, label='Direct Solar', ls='--')
        ax0.plot(dn_sol_flx,altitude,c='red', lw=2.0, label='Down Solar', ls='--')
        ax0.plot(up_sol_flx,altitude,c='grey', lw=2.0, label='Up Solar', ls='--')
        ax0.plot(up_th_flx,altitude,c='blue', lw=2.0, label='Up Thermal')
        ax0.plot(dn_th_flx,altitude,c='red', lw=2.0, label='Down Thermal')
        ax0.plot(0, 0, color='k', lw=2.0, label='Temperature')
        #ax0.axhspan(10,15, alpha=0.2, color='grey', label='Cloud')
        ax1 = ax0.twiny()
        ax1.plot(temperature, altitude, color='k', lw=2.0, label='Temperature [K]')
        ax1.set_xlabel(r'T [K]')
        ax0.set_xlabel(r"Flux [W m$^{-2}$]")
        ax0.set_ylabel("Altitude [km]")
        leg = ax0.legend()
        leg.get_frame().set_alpha(0.7)
        fig.show()
        if save:
            savetag = 'plot_hrt_'+ptitle+'.pdf'
            fig.savefig(savetag)
            print 'Saved plot: '+savetag

#################################################################################

def jacobians(path,jext,stream = 1, plot=False, getdata=True, xran=[1,5], pran=[1.e0,1.e5], ptitle='', radstep=False, contours=False, save=False, stitle=''):	

    '''Reads-in and parses SMART jacobians, returns a set of numpy arrays.
     
    Parameters
    ----------
    path : string
        SMART radiance output file
        Ex: 'earth_standard_hitran08_50_100000cm-1_clearsky_toa.rad'
    jext : string
        A unique indentifier used as the file extension for each Jacobian output --
        Specifies a physical quantity for which Jacobians have been calculated.
        Ex: '.j_O3' for Ozone Jacobians
    stream : int
        Specifices the upwelling stream to read-in

    Returns
    ----------
    jwl : float64
        Wavelength grid
    jwno : float64
        Wavenumber grid
    presprofile: float64
        Pressure grid
    jrad : float64
        Upwelling radiance grid 
    jacobians : float64
        Jacobian grid
    vprofile : float64
        Profile (or single quantity) corresponding to unperturbed state
    jfraction : float64
        Fractional change for which Jacobians were computed 

    Revision History
    ----------
    Written by J. Lustig-Yaeger June, 2015

    Examples    
    ----------  
    In [1]: from read_jacobians import get_jacobians
    In [2]: jext = '.j_O2'
    In [3]: dpath = '/astro/users/jlustigy/Models/smart_training/case06/'
    In [4]: radf = 'earth_standard_clearsky_jacobians_0.01_toa.rad'
    In [5]: jpath = dpath+radf
    In [6]: jwl, jwno, jpres, jrad, jacobians, vprofile, jfrac = get_jacobians(jpath,jext,stream=1)

    '''

    jfiles = ['_rad001_toa','_rad002_toa','_rad003_toa','_rad004_toa']
    jpath1 = path+jfiles[0]+jext
    jpath2 = path+jfiles[1]+jext
    jpath3 = path+jfiles[2]+jext
    jpath4 = path+jfiles[3]+jext
    jpaths = [jpath1,jpath2,jpath3,jpath4]
    
    if jext == '.j_pres' or jext == '.j_surf':
        single_stream = single_stream2
    else:
        single_stream = single_stream1
    
    if stream != 1 and stream != 2 and stream != 3 and stream != 4 and stream !='all':
        print('Invalid Stream. Must be 1, 2, 3, 4, or "all"')
        return
    elif stream == 1 or stream == 2 or stream == 3 or stream == 4:
        jwl, jwno, presprofile, jrad, jacobians, vprofile, jfraction \
        = single_stream(jpaths[stream-1])
    elif stream == 'all':
        jwl1, jwno1, presprofile1, jrad1, jacobians1, vprofile1, jfraction1 = \
        single_stream(jpath1)
        jwl2, jwno2, presprofile2, jrad2, jacobians2, vprofile2, jfraction2 = \
        single_stream(jpath2)
        jwl3, jwno3, presprofile3, jrad3, jacobians3, vprofile3, jfraction3 = \
        single_stream(jpath3)
        jwl4, jwno4, presprofile4, jrad4, jacobians4, vprofile4, jfraction4 = \
        single_stream(jpath4)
        
        jwl = np.array([jwl1,jwl2,jwl3,jwl4])
        jwno = np.array([jwno1,jwno2,jwno3,jwno4])
        presprofile = np.array([presprofile1,presprofile2,presprofile3,presprofile4])
        jrad = np.array([jrad1,jrad2,jrad3,jrad4])
        jacobians = np.array([jacobians1,jacobians2,jacobians3,jacobians4])
        vprofile = np.array([vprofile1,vprofile2,vprofile3,vprofile4])
        jfrac = jfraction
    
    if plot:
        wl, p, rad, jac, vp, jfrac = jwl, presprofile, jrad, jacobians, vprofile, jfraction
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        from matplotlib import gridspec
        from matplotlib import rc
        mpl.rcParams['font.size'] = 15.0    
        len1 = len(wl)
        mask = (wl > xran[0]) & (wl < xran[1])
        pmask = (p > pran[0]) & (p < pran[1])
        wl = wl[mask]
        rad = rad[mask]
        jac = jac[mask,:]
        jac = jac[:,pmask]
        p = p[pmask]
        vp = vp[pmask]
        len2 = len(wl)
        texstr = tex_molecule(jext[3:])
        if texstr == 'temp': texstr = 'T'
        if 'aer' in texstr: texstr = r'$\tau_{aero}$'
        percent = (len2*1.0/len1*1.0)*100.0
        print('Plotting '+str(len2)+' of '+str(len1)+' wavelength points ('+str(int(percent))+'%)')
    
        fig = plt.figure(figsize=(15,10))
        gs = gridspec.GridSpec(3, 2, height_ratios=[0.3,5, 2], width_ratios=[5,2]) 
        
        ax0 = plt.subplot(gs[1,0])
        vmag = np.max(np.abs(jac))
        pl = ax0.pcolor(wl, p, jac.T, cmap='seismic', vmin=-vmag, vmax=vmag)
        plt.title('')
        ax0.set_xlabel(r"Wavelength [$\mu$m]")
        ax0.set_ylabel("Pressure [Pa]")
        ax0.invert_yaxis()
        ax0.semilogy()
        ax0.set_xlim(left=xran[0],right=xran[1])
        ax0.axes.get_xaxis().set_visible(False)

        if contours:
            cnt1 = ax0.contour(wl, p, jac.T, colors='k')
            ax0.clabel(cnt1, inline=3, fontsize=14,fmt='%1.1f', colors='k')

        # Add color bar using gridspec
        cbaxes = plt.subplot(gs[0,0]) 
        cb = plt.colorbar(pl, cax = cbaxes, orientation="horizontal", use_gridspec=True)
        cb.set_label(texstr+r" Jacobians [dR/d"+texstr+"]$ _{ \Delta = "+str(jfrac)+" } $", labelpad=-90)
        cbaxes.xaxis.tick_top()
        
        # attempting to modify radiance with jacobians
        newrad = rad + np.inner(jac,vp)
        vp_plus = vp + vp*jfrac
        vp_minus = vp - vp*jfrac
        
        ax1 = plt.subplot(gs[2,0],sharex=ax0)
        ax1.plot(wl, rad, color='k')
        if radstep:
            ax1.plot(wl, newrad, color='k', ls='--')
        ax1.plot(wl, rad, color='k')
        ax1.set_xlim(left=xran[0],right=xran[1])
        ax1.set_ylabel(r"Radiance") #[W m$^{-2}$ $\mu$m$^{-1}$ sr$^{-1}$]
        ax1.set_xlabel(r"Wavelength [$\mu$m]")
        #ax1.semilogx()
        #ax2 = ax1.twinx()
        #ax2.plot(wl2,qe)
    
        ax2 = plt.subplot(gs[1,1])
        ax2.plot(vp, p*1e-5, color='k')
        if radstep:
            ax2.fill_betweenx(p, vp_minus, vp_plus, color='k', alpha='0.5')
        #ax2.set_xlim(left=wlrange[0],right=wlrange[1])
        ax2.set_ylabel(r"Pressure [bars]")
        ax2.set_xlabel(r"Mixing Ratio")
        if jext == '.j_temp':
            ax2.set_xlabel(r"T [K]")
        elif 'aer' in jext:
            ax2.set_xlabel(r"$\tau_{aero}$")
            ax2.semilogx()
        else:
            ax2.set_xlabel(texstr+" Mixing Ratio")
            ax2.semilogx()
        #ax2.axes.get_yaxis().set_visible(False)
        ax2.invert_yaxis()
        ax2.yaxis.tick_right() 
        ax2.yaxis.set_label_position("right")
        #ax2.set_xlim(left=-0.000001, right=np.amax(vp)+0.0001)
        #ax2.set_xticks(np.arange(-0.01,np.amax(vp)+0.01))
        ax2.semilogy()
        #ax2 = ax1.twinx()
        #ax2.plot(wl2,qe)

        fig.subplots_adjust(hspace=0.1)
        fig.subplots_adjust(wspace=0.05)

        if save:
            savetag = 'plot_jacobians_'+stitle+'_'+ptitle+'.pdf'
            fig.savefig(savetag)
            print 'Saved plot: '+savetag
    
        fig.show();

        
    if getdata:
        return jwl, jwno, presprofile, jrad, jacobians, vprofile, jfraction
    else:
        return

    
def read_header(arr):
    arr[0:30]
    levels = arr[0][0]
    vlength = arr[0][2]
    jfraction = arr[1][0]
    vprofile = arr[1][1:]
    i=2
    while len(vprofile) < vlength:
        vprofile = np.concatenate([vprofile,arr[i]])
        i += 1
    wnorange = [arr[i][0],arr[i][1]]
    presprofile = arr[i][2:]
    i += 1
    while len(presprofile) < levels:
        presprofile = np.concatenate([presprofile,arr[i]])
        i += 1
    return i, levels, jfraction, vprofile, wnorange, presprofile

def read_header2(arr):
    layers = arr[0][0]
    jfraction = arr[1][0]
    jquant = arr[1][1]
    i=2

    wnorange = [arr[i][0],arr[i][1]]
    presprofile = arr[i][2:]
    i += 1
    while len(presprofile) < layers:
        presprofile = np.concatenate([presprofile,arr[i]])
        i += 1
    return i, layers, jfraction, jquant, wnorange, presprofile

def parse_jacobians(arr):
    
    # Extract info from header
    index, layers, jfraction, vprofile, wnorange, presprofile \
    = read_header(arr)
    
    # Define a bunch of quantities from header of file
    x = len(arr)*1.0
    y = len(arr[index])*1.0
    p_tot = len(vprofile)*1.0
    z = np.mod(p_tot+3, y)*1.0
    a = np.floor((p_tot+3)/y)
    
    # if there is an additional unfilled row then d=1, else d=0
    if z == 0: 
        d = 0*1.0
    else:
        d = 1.0
    
    h = int(a+d) # rows in file per true row
    N = (x-index)/(a+d) # total number of true rows in file
    wno1 = arr[index][1]*1.0
    wno2 = arr[index+int(a+d)][1]*1.0
    ranwno = (wnorange[1] - wnorange[0])*1.0
    deltawno = (wno2 - wno1)*1.0
    wnotot = ranwno / deltawno
    
    # Construct structured array for data (N, layers+3)
    J = np.zeros((N,p_tot+3))
    
    # loop over all N in J, filling true row with appropriate file rows
    for i in range(int(N)):
        try:
            J[i,:] = np.hstack(arr[(h*i)+index:(h*(i+1)+index)])
        except IndexError:
            print("IndexError: i = "+str(i), J.shape) 
            break
    
    # decompose J into discrete pieces
    jwl, jwno, jrad, jacobians = J[:,0], J[:,1], J[:,2], J[:,3:]
    
    return jwl, jwno, presprofile, jrad, jacobians, vprofile, jfraction

def parse_jacobians2(arr):
    
    # Extract info from header
    index, layers, jfraction, jquant, wnorange, presprofile \
    = read_header2(arr)
        
    # Construct structured array for data 
    J = np.zeros((len(arr[index:]),4))
    Jarr= np.array(arr[index:])
    
    # loop over all indexes assigning each row to J
    for i in range(len(Jarr)):
        J[i,0] = Jarr[i][0]
        J[i,1] = Jarr[i][1]
        J[i,2] = Jarr[i][2]
        J[i,3] = Jarr[i][3]

    jwl, jwno, jrad, jacobians = J[:,0], J[:,1], J[:,2], J[:,3]
    
    return jwl, jwno, presprofile, jrad, jacobians, jquant, jfraction

def single_stream1(path):
    
    #print("Opening j_file: "+path)
    arrays = [np.array(map(float, line.split())) for line in open(path)]
    arrays = np.array(arrays)
    
    #print("Length of file: "+str(len(arrays)))
    
    jwl1, jwno1, presprofile1, jrad1, jacobians1, vprofile1, jfraction1 = \
    parse_jacobians(arrays)
    
    return jwl1, jwno1, presprofile1, jrad1, jacobians1, vprofile1, jfraction1

def single_stream2(path):
    
    #print("Opening j_file: "+path)
    arrays = [np.array(map(float, line.split())) for line in open(path)]
    arrays = np.array(arrays)
    
    #print("Length of file: "+str(len(arrays)))
    
    jwl1, jwno1, presprofile1, jrad1, jacobians1, jquant1, jfraction1 = \
    parse_jacobians2(arrays)
    
    return jwl1, jwno1, presprofile1, jrad1, jacobians1, jquant1, jfraction1

#################################################################################

def plot_jcolor(wl, p, rad, jac, vp, jfrac=.01, ptitle='', wlrange=[1.0,2.0], radstep=False, contours=False,save=False, stitle=''):

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib import gridspec
    from matplotlib import rc
    mpl.rcParams['font.size'] = 20.0    

    len1 = len(wl)
    
    mask = (wl > wlrange[0]) & (wl < wlrange[1])
    wl = wl[mask]
    rad = rad[mask]
    jac = jac[mask,:]
    
    len2 = len(wl)
    percent = (len2*1.0/len1*1.0)*100.0
    print('Plotting '+str(len2)+' of '+str(len1)+' wavelength points ('+str(int(percent))+'%)')
    
    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(3, 2, height_ratios=[0.3,5, 2], width_ratios=[5,2]) 
    ax0 = plt.subplot(gs[1,0])
    pl = ax0.pcolor(wl, p, jac.T)
    plt.title('')
    ax0.set_xlabel(r"Wavelength [$\mu$m]")
    ax0.set_ylabel("Pressure [Pa]")
    ax0.invert_yaxis()
    ax0.semilogy()
    ax0.set_xlim(left=wlrange[0],right=wlrange[1])
    ax0.axes.get_xaxis().set_visible(False)

    if contours:
        cnt1 = ax0.contour(wl, p, jac.T, colors='k')
        ax0.clabel(cnt1, inline=3, fontsize=14,fmt='%1.1f', colors='k')

    # Add color bar using gridspec
    cbaxes = plt.subplot(gs[0,0]) 
    cb = plt.colorbar(pl, cax = cbaxes, orientation="horizontal", use_gridspec=True)
    cb.set_label(ptitle+r" Jacobians [$dR / d\mathbf{x}_i$]", labelpad=-70)
    cbaxes.xaxis.tick_top()
    
    # attempting to modify radiance with jacobians
    newrad = rad + np.inner(jac,vp)
    vp_plus = vp + vp*jfrac
    vp_minus = vp - vp*jfrac
    
    ax1 = plt.subplot(gs[2,0],sharex=ax0)
    ax1.plot(wl, rad, color='k')
    if radstep:
        ax1.plot(wl, newrad, color='k', ls='--')
    ax1.plot(wl, rad, color='k')
    ax1.set_xlim(left=wlrange[0],right=wlrange[1])
    ax1.set_ylabel(r"Radiance") #[W m$^{-2}$ $\mu$m$^{-1}$ sr$^{-1}$]
    ax1.set_xlabel(r"Wavelength [$\mu$m]")
    #ax1.semilogx()
    #ax2 = ax1.twinx()
    #ax2.plot(wl2,qe)
    
    ax2 = plt.subplot(gs[1,1],sharey=ax0)
    ax2.plot(vp, p, color='k')
    if radstep:
        ax2.fill_betweenx(p, vp_minus, vp_plus, color='k', alpha='0.5')
    #ax2.set_xlim(left=wlrange[0],right=wlrange[1])
    ax2.set_ylabel(r"Pressure [Pa]")
    ax2.set_xlabel(r"Mixing Ratio")
    if ptitle == 'Temperature':
        ax2.set_xlabel(r"T [K]")
    else:
         ax2.semilogx()
    ax2.axes.get_yaxis().set_visible(False)
    #ax2.set_xlim(left=-0.000001, right=np.amax(vp)+0.0001)
    #ax2.set_xticks(np.arange(-0.01,np.amax(vp)+0.01))

    #ax2.semilogx()
    #ax2 = ax1.twinx()
    #ax2.plot(wl2,qe)

    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.05)

    if save:
        fig.savefig('plot_jcolor_'+stitle+'_'+ptitle+'.pdf')
    
    fig.show();

    return

def tex_molecule(formula):
    new = r''
    for char in formula:
        if char.isdigit():
            tmp = '$_{'+char.upper()+'}$'
        else:
            tmp = char
        new = new+tmp
    return new




