import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

class Filter(object):
    """Filter for telescope imaging mode.
    
    Parameters
    ----------
    name : string
        Name of filter
    bandcenter : float
        Wavelength at bandcenter (um)
    FWHM : float
        Full-width as half-maximum (um)
    wl : array
        Wavelength grid for response function (um)
    response : array
        Filter response function
    notes : string
        Notes of filter, e.g. 'Johnson-Cousins'
    """
    
    def __init__(self, name=None, bandcenter=None, FWHM=None, wl=None, response=None, notes=''):
        self.name=name
        self.bandcenter=bandcenter
        self.FWHM=FWHM
        self.wl=wl
        self.response=response
        self.notes=notes
        
    def __str__(self):
        string = 'Filter: \n------------\n'+\
            '- Name : '+"%s" % (self.name)+'\n'+\
            '- Bandcenter (um)  : '+"%s" % (self.bandcenter)+'\n'+\
            '- Full width at half max (um)  : '+"%s" % (self.FWHM)+'\n'+\
            '- Wavelength array (um)  : '+"%s" % (self.wl)+'\n'+\
            '- Response : '+"%s" % (self.response)+'\n'+\
            '- Notes  : '+"%s" % (self.notes)+'\n'
        return string

class Wheel(object): 
    """Filter Wheel. Contains different filters as attributes.
    """
    def __init__(self):
        pass
    
    def add_new_filter(self, filt, name='new_filter'):
        """Adds new filter to wheel
        
        Parameters
        ----------
        filt : Filter
            New Filter object to be added to wheel
        name : string (optional)
            Name to give new filter attribute
        """
        setattr(self, name, filt)
     
    def plot(self, ax=None):
        
        if ax == None:
            fig = plt.figure(figsize=(14,10))
            gs = gridspec.GridSpec(1,1) 
            ax1 = plt.subplot(gs[0])
            ax1.set_ylabel(r"Filter Response")
            ax1.set_xlabel(r"Wavelength [$\mu$m]")
            #ax1.set_ylim([0.0,1.0])
        else:
            ax1 = ax.twinx()
            #ax1.set_ylim([0.0,10.0])
            ax1.axes.get_yaxis().set_visible(False)          

        Nfilt = len(self.__dict__)   
        #colors,scalarMap,cNorm = scalarmap(np.arange(Nfilt),cmap='Dark2')     
        i = 0
        fmax = 1.0
        for attr, value in self.__dict__.iteritems():
            if np.max(value.response) > fmax: fmax = np.max(value.response)
            wl, response = value.wl, value.response
            #ax1.plot(wl,response, lw=3.0, label=value.name, c=colors[i])
            ax1.fill_between(wl,response, color='purple', alpha=0.3)
            i += 1
        ax1.set_ylim([0.0,fmax*10.])
        
        if ax==None:
            return ax1
        
    def __str__(self):
        string = []
        tdict = sorted(self.__dict__.iteritems(), key=lambda x: x[1].bandcenter)
        for attr, value in tdict:
            string.append(attr)
        print string 
        return ''#str(string)    


def read_jc():
    path = 'UBVRI/'
    # set file path relative to this file
    path = os.path.join(os.path.dirname(__file__), path)
    U = np.genfromtxt(path+'bess-u.pass')
    U[:,0] = U[:,0]/1e4
    B = np.genfromtxt(path+'bess-b.pass')
    B[:,0] = B[:,0]/1e4
    V = np.genfromtxt(path+'bess-v.pass')
    V[:,0] = V[:,0]/1e4
    R = np.genfromtxt(path+'bess-r.pass')
    R[:,0] = R[:,0]/1e4
    I = np.genfromtxt(path+'bess-i.pass')
    I[:,0] = I[:,0]/1e4
    filters = np.array([U,B,V,R,I])
    filter_names = ['U','B','V','R','I']
    bandcenters = np.array([365.6, 435.3, 547.7, 634.9, 800.0]) / 1e3
    FWHM = np.array([34.0, 78.1, 99.1, 106.56, 289.2]) / 1e3
    return filters, filter_names, bandcenters, FWHM

class johnson_cousins(Wheel):
    
    def __init__(self):
        
        filters, filter_names, bandcenters, FWHM = read_jc()
        
        self.U=Filter(name='U', bandcenter=bandcenters[0], FWHM=FWHM[0], wl=filters[0][:,0], response=filters[0][:,1], notes='Johnson-Cousins')
        self.B=Filter(name='B', bandcenter=bandcenters[1], FWHM=FWHM[1], wl=filters[1][:,0], response=filters[1][:,1], notes='Johnson-Cousins')
        self.V=Filter(name='V', bandcenter=bandcenters[2], FWHM=FWHM[2], wl=filters[2][:,0], response=filters[2][:,1], notes='Johnson-Cousins')
        self.R=Filter(name='R', bandcenter=bandcenters[3], FWHM=FWHM[3], wl=filters[3][:,0], response=filters[3][:,1], notes='Johnson-Cousins')
        self.I=Filter(name='I', bandcenter=bandcenters[4], FWHM=FWHM[4], wl=filters[4][:,0], response=filters[4][:,1], notes='Johnson-Cousins')

        
def read_landsat():
    path = 'LANDSAT/'
    # set file path relative to this file
    path = os.path.join(os.path.dirname(__file__), path)
    blue = np.genfromtxt(path+'Blue.txt', skip_header=1)
    green = np.genfromtxt(path+'Green.txt', skip_header=1)
    red = np.genfromtxt(path+'Red.txt', skip_header=1)
    coastal = np.genfromtxt(path+'CostalAerosol.txt', skip_header=1)
    cirrus = np.genfromtxt(path+'Cirrus.txt', skip_header=1)
    nir = np.genfromtxt(path+'NIR.txt', skip_header=1)
    pan = np.genfromtxt(path+'Pan.txt', skip_header=1)
    swir1 = np.genfromtxt(path+'SWIR1.txt', skip_header=1)
    swir2 = np.genfromtxt(path+'SWIR2.txt', skip_header=1)
    LANDSAT_names = ['Coastal Aerosols','Blue','Green','Red','NIR','SWIR1','SWIR2','Pan','Cirrus']
    titles = ['wl','response','std','bandwidth','FWHM_low','FWHM_high','bandcenter']
    wl = [coastal[:,0]/1e3, blue[:,0]/1e3, green[:,0]/1e3, red[:,0]/1e3, nir[:,0]/1e3, swir1[:,0]/1e3, swir2[:,0]/1e3, pan[:,0]/1e3, cirrus[:,0]/1e3]
    response = [coastal[:,1], blue[:,1], green[:,1], red[:,1], nir[:,1], swir1[:,1], swir2[:,1], pan[:,1], cirrus[:,1]]
    FWHM = np.array([15.98, 60.04, 57.33, 37.47, 28.25, 84.72, 186.66, 172.40, 20.39]) / 1e3
    bandcenters = np.array([442.96, 482.04, 561.41, 654.59, 864.67, 1608.86, 2200.73, 589.50, 1373.43]) / 1e3
    return wl, response, LANDSAT_names, FWHM, bandcenters

class landsat(Wheel):
    
    def __init__(self):
        
        wl, response, LANDSAT_names, FWHM, bandcenters = read_landsat()
        
        self.CA=Filter(name=LANDSAT_names[0], bandcenter=bandcenters[0], FWHM=FWHM[0], wl=wl[0], response=response[0], notes='LANDSAT')
        self.B=Filter(name=LANDSAT_names[1], bandcenter=bandcenters[1], FWHM=FWHM[1], wl=wl[1], response=response[1], notes='LANDSAT')
        self.G=Filter(name=LANDSAT_names[2], bandcenter=bandcenters[2], FWHM=FWHM[2], wl=wl[2], response=response[2], notes='LANDSAT')
        self.R=Filter(name=LANDSAT_names[3], bandcenter=bandcenters[3], FWHM=FWHM[3], wl=wl[3], response=response[3], notes='LANDSAT')
        self.NIR=Filter(name=LANDSAT_names[4], bandcenter=bandcenters[4], FWHM=FWHM[4], wl=wl[4], response=response[4], notes='LANDSAT')
        self.SWIR1=Filter(name=LANDSAT_names[5], bandcenter=bandcenters[5], FWHM=FWHM[5], wl=wl[5], response=response[5], notes='LANDSAT')
        self.SWIR2=Filter(name=LANDSAT_names[6], bandcenter=bandcenters[6], FWHM=FWHM[6], wl=wl[6], response=response[6], notes='LANDSAT')
        self.Pan=Filter(name=LANDSAT_names[7], bandcenter=bandcenters[7], FWHM=FWHM[7], wl=wl[7], response=response[7], notes='LANDSAT')
        self.Cirrus=Filter(name=LANDSAT_names[8], bandcenter=bandcenters[8], FWHM=FWHM[8], wl=wl[8], response=response[8], notes='LANDSAT')
        
        
def read_jc2():
    path = 'UBVRI2/'
    # set file path relative to this file
    path = os.path.join(os.path.dirname(__file__), path)
    U = np.genfromtxt(path+'Bessel_U-1.txt', skip_header=1)
    U[:,0] = U[:,0]/1e3
    B = np.genfromtxt(path+'Bessel_B-1.txt', skip_header=1)
    B[:,0] = B[:,0]/1e3
    V = np.genfromtxt(path+'Bessel_V-1.txt', skip_header=1)
    V[:,0] = V[:,0]/1e3
    R = np.genfromtxt(path+'Bessel_R-1.txt', skip_header=1)
    R[:,0] = R[:,0]/1e3
    I = np.genfromtxt(path+'Bessel_I-1.txt', skip_header=1)
    I[:,0] = I[:,0]/1e3
    filters = np.array([U[::-1,:],B[::-1,:],V[::-1,:],R[::-1,:],I[::-1,:]])
    filter_names = ['U','B','V','R','I']
    bandcenters = np.array([365.6, 435.3, 547.7, 634.9, 879.7]) / 1e3
    FWHM = np.array([34.0, 78.1, 99.1, 106.56, 289.2]) / 1e3
    return filters, filter_names, bandcenters, FWHM

class johnson_cousins2(Wheel):
    
    def __init__(self):
        
        filters, filter_names, bandcenters, FWHM = read_jc2()
        
        self.U=Filter(name='U', bandcenter=bandcenters[0], FWHM=FWHM[0], wl=filters[0][:,0], response=filters[0][:,1], notes='Johnson-Cousins')
        self.B=Filter(name='B', bandcenter=bandcenters[1], FWHM=FWHM[1], wl=filters[1][:,0], response=filters[1][:,1], notes='Johnson-Cousins')
        self.V=Filter(name='V', bandcenter=bandcenters[2], FWHM=FWHM[2], wl=filters[2][:,0], response=filters[2][:,1], notes='Johnson-Cousins')
        self.R=Filter(name='R', bandcenter=bandcenters[3], FWHM=FWHM[3], wl=filters[3][:,0], response=filters[3][:,1], notes='Johnson-Cousins')
        self.I=Filter(name='I', bandcenter=bandcenters[4], FWHM=FWHM[4], wl=filters[4][:,0], response=filters[4][:,1], notes='Johnson-Cousins')
