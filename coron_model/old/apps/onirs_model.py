#####BOKEH O/NIRS NOISE MODEL SIMULATOR#####
# This code produces an interactive browser widget that runs
# the o/nirs noise model
# This was written based on the coronagraph noise model code.
# Warning: Some defunct variables may still be present
# from the coronagraph noise model.
#
# To run this code on your local machine, type
# bokeh serve --show onirs_model.py
# 
################################################

# Import some standard python packages

import numpy as np
from astropy.io import fits, ascii 
import pdb
import sys
import os 
from astropy.table import Table, Column
import os
from bokeh.io import curdoc
from bokeh.client import push_session
import pysynphot as S

from bokeh.themes import Theme 
import yaml 
from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, HBox, VBoxForm, HoverTool, Paragraph, Range1d, DataRange1d, Label, DataSource
from bokeh.models.glyphs import Text
from bokeh.layouts import column, row, WidgetBox 
from bokeh.models.widgets import Slider, Panel, Tabs, Div, TextInput, RadioButtonGroup, Select, RadioButtonGroup
from bokeh.io import hplot, vplot, curdoc, output_file, show, vform
from bokeh.models.callbacks import CustomJS
from bokeh.embed import components, autoload_server


import coronagraph as cg  # Import coronagraph (O/NIRS)  model

#allow it to run it from other folders and still know where planet folder is
planetdir = "../coronagraph/planets/" #new path compared to before
relpath = os.path.join(os.path.dirname(__file__), planetdir)

################################
# PARAMETERS
################################

# Integration time (hours)
Dt = 24.0 # - SLIDER

# Telescopes params
diam = 10. # mirror diameter - SLIDER
Res = 100. # vis resolution - SLIDER
Tsys = 270. # system temperature - SLIDER

#Object parameters
#Ro = 1.         # object radius (solar radii) - SLIDER
#d = 1./206265. # 1 AU
#d    = 10.      # no - we are going to have the user specify magnitude because distance is a ridiculous quantity when it ranges from solar system-like to extragalactic scales

mag = 10. #nominal magnitude for sun spectrum

# Instrumental Params
owa = 3000. #OWA scaling factor - SLIDER
iwa = 0.001 #IWA scaling factor - SLIDER
De_UV = 1e-4  # dark current - not making slider for now
De_VIS = 1e-4 # dark current - not making slider for now
De_NIR = 1e-4 # dark current - not making slider for now
Re_UV = 0.01   # read noise - not making slider for now 
Re_VIS = 0.01  # read noise - not making slider for now 
Re_NIR = 0.01  # read noise - not making slider for now 
Res = 100
Res_UV = 20
Res_NIR = 150
Re = 0.1 # read noise - 
Dtmax = 1.0 # max single exposure time - SLIDER
wantsnr = 10. #for exposure time calculator - SLIDER

#unit conversion things
erg_s_per_watt = 1e7 
cm2_per_m2 = 1e4
A_per_micron = 1e4 


# Template
template = ''
global template
global comparison
global Teff
global Ts
global erg_s_per_watt
global cm2_per_m2
global A_per_micron
global De_UV
global De_VIS
global De_NIR
global Re_UV
global Re_VIS
global Re_NIR

################################
# READ-IN DATA

# Read-in Sun spectrum file to start 

fn = 'earth_quadrature_radiance_refl.dat'
fn = os.path.join(relpath, fn)
model = np.loadtxt(fn, skiprows=8)
lamhr = model[:,0]
radhr = model[:,1]
Fohr = model[:,2] #* (1.495e11)**2 / (6.95e8)**2 # convert to flux @ stellar surface
#print 'Fohr =', Fohr
# Calculate hi-resolution reflectivity
Fohr_bb = cg.noise_routines.Fstar(lamhr, 5777, 1., 1., AU=True) # stellar flux (comparison)
lammin = 0.4
lammax = 5.
planet_label = ['']
#pdb.set_trace()
Fohr_ = Fohr
lamhr_ = lamhr

# read in Vega spectrum for flux calibration 
sp = S.Vega
VegaWave = S.Vega.wave #angstroms
VegaFlux = S.Vega.flux

#convert to the input magnitude
# erg/s/cm2/A -> W/m2/um
VegaFlux = VegaFlux / erg_s_per_watt * cm2_per_m2 * A_per_micron
VegaWave = VegaWave / A_per_micron
global VegaFlux
global VegaWave


#V band
Vmin = 0.5 #microns
Vmax = 0.58 #microns
global Vmin
global Vmax


#integrate the visible flux
vegaindex = np.where((VegaWave > Vmin) & (VegaWave < Vmax))
objindex = np.where((lamhr > Vmin) & (lamhr < Vmax))
VegaIntFlux = np.trapz(VegaFlux[vegaindex], x=VegaWave[vegaindex])
ObjIntFlux = np.trapz(Fohr_[objindex], x=lamhr[objindex])
F_mag = 10**(-0.4*(mag))*VegaIntFlux #the flux it WOULD have at the user-specified magnitude
factor = F_mag / ObjIntFlux
Fohr_ = Fohr_ * factor
#import pdb; pdb.set_trace()



################################
# RUN O/NIRS MODEL
################################

# Run onirs with default LUVOIR telescope 
lam, dlam, Fo, q, co, cz, cD, cR, cth, DtSNR = \
    cg.count_rates_onirs(Fohr_, lamhr, lammin=lammin, lammax=lammax, Res=Res, Res_UV=Res_UV, Res_NIR=Res_NIR, diam=diam, Tsys=Tsys,De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, Dtmax=Dtmax, GROUND=False, THERMAL=True,  wantsnr=wantsnr)

print 'co =', co

# Calculate background photon count rates
cb = (cz + cD + cR + cth)
# Convert hours to seconds
Dts = Dt * 3600.
# Calculate signal-to-noise assuming background subtraction (the "2")
SNR  = co*Dts/np.sqrt((co + 2*cb)*Dts)
# Calculate 1-sigma errors
sig= Fo/SNR
# Add gaussian noise to flux ratio
spec = Fo + np.random.randn(len(Fo))*sig #want to plot planet flux, not counts

#update params
lastlam = lam
snr_ymax = np.max(Fo)
yrange=[snr_ymax]


#data
theobject = ColumnDataSource(data=dict(lam=lam, Fo=Fo, spec=spec, downerr=(spec-sig), uperr=(spec+sig), cz=cz*Dts,  cD=cD*Dts, cR=cR*Dts, cth=cth*Dts, co=co*Dts))
expobject = ColumnDataSource(data=dict(lam=lam[np.isfinite(DtSNR)], DtSNR=DtSNR[np.isfinite(DtSNR)])) 
plotyrange = ColumnDataSource(data = dict(yrange=yrange))
textlabel = ColumnDataSource(data=dict(label = planet_label))


################################
# BOKEH PLOTTING
################################
#plots spectrum and exposure time
snr_plot = Figure(plot_height=500, plot_width=750, 
                  tools="crosshair,pan,reset,resize,save,box_zoom,wheel_zoom,hover",
                  toolbar_location='right', x_range=[0.4, 5.0], y_range=[min(Fo)*0.9, max(Fo)*1.1])

snr_plot.background_fill_color = "beige"
snr_plot.background_fill_alpha = 0.5
snr_plot.yaxis.axis_label='Flux [W/m**2/s]' 
snr_plot.xaxis.axis_label='Wavelength [micron]'
snr_plot.title.text = 'Object Spectrum: Sun' #initial spectrum is Sun at 10 pc

snr_plot.line('lam','Fo',source=theobject,line_width=2.0, color="darkgreen", alpha=0.7)
snr_plot.circle('lam', 'spec', source=theobject, fill_color='lightgreen', line_color='black', size=8) 
snr_plot.segment('lam', 'downerr', 'lam', 'uperr', source=theobject, line_width=1, line_color='grey', line_alpha=0.5)


#text on plot
glyph = Text(x=0.25, y=-0.19, text="label", text_font_size='9pt', text_font_style='bold', text_color='blue')
#attempting to outline the text here for ease of visibility... 
glyph2 = Text(x=0.245, y=-0.19, text="label", text_font_size='9pt', text_font_style='bold', text_color='white')
glyph3 = Text(x=0.25, y=-0.195, text="label", text_font_size='9pt', text_font_style='bold', text_color='white')
glyph4 = Text(x=0.25, y=-0.845, text="label", text_font_size='9pt', text_font_style='bold', text_color='white')
glyph5 = Text(x=0.255, y=-0.19, text="label", text_font_size='9pt', text_font_style='bold', text_color='white')
snr_plot.add_glyph(textlabel, glyph2)
snr_plot.add_glyph(textlabel, glyph3)
snr_plot.add_glyph(textlabel, glyph4)
snr_plot.add_glyph(textlabel, glyph5)
snr_plot.add_glyph(textlabel, glyph)

#hovertool
hover = snr_plot.select(dict(type=HoverTool))
from math import log
#logobj = [log(y,10) for y in co]
#logobj = [round(y,4) for y in logobj]
#logobj = np.asarray(logobj)
hover.tooltips = [
   ('object', '@co{int}'),
   ('zodi', '@cz{int}'),
   ('dark current', '@cD{int}'),
   ('read noise', '@cR{int}'),
   ('thermal', '@cth{int}')
]

ptab1 = Panel(child=snr_plot, title='Spectrum')
#ptab2 = Panel(child=exp_plot, title='Exposure Time')
ptabs = Tabs(tabs=[ptab1])
show(ptabs)

################################
#  PROGRAMS
################################

def change_filename(attrname, old, new): 
   format_button_group.active = None 


instruction0 = Div(text="""Specify a filename here:
                           (no special characters):""", width=300, height=15)
text_input = TextInput(value="filename", title=" ", width=100)
instruction1 = Div(text="""Then choose a file format here:""", width=300, height=15)
format_button_group = RadioButtonGroup(labels=["txt", "fits"])
instruction2 = Div(text="""The link to download your file will appear here:""", width=300, height=15)
link_box  = Div(text=""" """, width=300, height=15)


def i_clicked_a_button(new): 
    filename=text_input.value + {0:'.txt', 1:'.fits'}[format_button_group.active]
    print "Your format is   ", format_button_group.active, {0:'txt', 1:'fits'}[format_button_group.active] 
    print "Your filename is: ", filename 
    fileformat={0:'txt', 1:'fits'}[format_button_group.active]
    link_box.text = """Working""" 
 
    t = Table(planet.data)
    t = t['lam', 'spec','cratio','uperr','downerr'] 

    if (format_button_group.active == 1): t.write(filename, overwrite=True) 
    if (format_button_group.active == 0): ascii.write(t, filename)
 
    os.system('gzip -f ' +filename) 
    os.system('cp -rp '+filename+'.gz /home/jtastro/jt-astro.science/outputs') 
    print    """Your file is <a href='http://jt-astro.science/outputs/"""+filename+""".gz'>"""+filename+""".gz</a>. """

    link_box.text = """Your file is <a href='http://jt-astro.science/outputs/"""+filename+""".gz'>"""+filename+""".gz</a>. """


#########################################
# GET DATA FROM USER AND UPDATE PLOT
#########################################

def update_data(attrname, old, new):
    print 'Updating model for exptime = ', exptime.value, ' with magnitude = ', magnitude.value 
    print 'diameter (m) = ', diameter.value, 'resolution = ', resolution.value
    print '                   temperature (K) = ', temperature.value
    print 'You have chosen object spectrum: ', template.value
    try:
       lasttemplate
    except NameError:
       lasttemplate = 'Sun' #default first spectrum
       
    global lasttemplate
    global lamhr_
    global Teff_
    global radhr_
    global Fohr_

    
# Read-in new spectrum file only if changed
    if template.value != lasttemplate:
       if template.value == 'Sun':
          #already in W/m2/um so don't need to convert
          fn = 'earth_quadrature_radiance_refl.dat'
          fn = os.path.join(relpath, fn)
          model = np.loadtxt(fn, skiprows=8)
          lamhr_ = model[:,0]
          radhr_ = model[:,1]
          Fohr_ = model[:,2]
#          print 'Fohr_ =', Fohr_
          lammin = 0.4
          lammax = 5.
          planet_label = ['']
          
          
       global lammin
       global lammax
       global planet_label
       lammin=min(lamhr_)
       if lammin <= 0.4:
          lammin = 0.4
       lammax=5.

    #convert to user-specified magnitude
    vegaindex = np.where((VegaWave > Vmin) & (VegaWave < Vmax))
    objindex = np.where((lamhr_ > Vmin) & (lamhr_ < Vmax))
    VegaIntFlux = np.trapz(VegaFlux[vegaindex], x=VegaWave[vegaindex])
    ObjIntFlux = np.trapz(Fohr_[objindex], x=lamhr_[objindex])
    F_mag = 10**(-0.4*(magnitude.value))*VegaIntFlux #the flux it WOULD have at the user-specified magnitude
    factor = F_mag / ObjIntFlux
    Fohr_ = Fohr_ * factor

    print "ground based = ", ground_based.value
    if ground_based.value == "No":
       ground_based_ = False
    if ground_based.value == "Yes":
       ground_based_ = True
    
    # Run o/nirs 
    lam, dlam, Fo, q, co, cz, cD, cR, cth, DtSNR = \
        cg.count_rates_onirs(Fohr_, lamhr_,  Res=resolution.value, Tsys=temperature.value, lammin=lammin, lammax=lammax,De_UV=De_UV, De_VIS=De_VIS, De_NIR=De_NIR, Re_UV=Re_UV, Re_VIS=Re_VIS, Re_NIR=Re_NIR, diam=diameter.value, Dtmax = dtmax.value, THERMAL=True, GROUND=ground_based_, wantsnr=want_snr.value, gain=gain.value, ntherm=ntherm.value)

 #   lam, dlam, Fo, q, co, cz, cD, cR, cth, DtSNR = \
  #      cg.count_rates_onirs(Fohr_, lamhr, Ro, d, lammin=lammin, lammax=lammax, Res=Res, diam=diam, Tsys=Tsys,De=De, Re=Re, Dtmax=Dtmax, GROUND=False, THERMAL=True,  wantsnr=wantsnr)

  
    # Calculate background photon count rates
    cb = (cz + cD + cR + cth)
    # Convert hours to seconds
    Dts = exptime.value * 3600.
    # Calculate signal-to-noise assuming background subtraction (the "2")
    SNR  = co*Dts/np.sqrt((co + 2*cb)*Dts)
    # Calculate 1-sigma errors
    sig= Fo/SNR
    # Add gaussian noise to flux ratio
    spec = Fo + np.random.randn(len(Fo))*sig #want to plot planet flux, not counts

    lastlam = lam
    lastFo = Fo
    global lastlam


    #UPDATE DATA
    #data
    theobject.data = dict(lam=lam, Fo=Fo, spec=spec, downerr=(spec-sig), uperr=(spec+sig), cz=cz*Dts,  cD=cD*Dts, cR=cR*Dts, cth=cth*Dts, co=co*Dts)
    plotyrange.data = dict(yrange=yrange)
    textlabel.data = dict(label = planet_label)

    
    #planet.data = dict(lam=lam, cratio=Cratio*1e9, spec=spec*1e9, downerr=(spec-sig)*1e9, uperr=(spec+sig)*1e9, cz=cz*Dts, cez=cez*Dts, csp=csp*Dts, cD=cD*Dts, cR=cR*Dts, cth=cth*Dts, cp=cp*Dts)
    #expplanet.data = dict(lam=lam[np.isfinite(DtSNR)], DtSNR=DtSNR[np.isfinite(DtSNR)])
     #make the data the time for a given SNR if user wants this:
   # textlabel.data = dict(label=planet_label)

    format_button_group.active = None
    lasttemplate = template.value

         
    #######PLOT UPDATES#######    
    global snr_ymax_
    global snr_ymin_

   # ii = np.where(lam < 2.5) #only want where reflected light, not thermal
  #  iii = np.where(lamC < 2.5)  #only want where reflected light, not thermal
   # pdb.set_trace()
    #Cratio_ok = Cratio[ii]
    #CratioC_ok = CratioC[iii]
    #Cratio_ok = Cratio_ok[~np.isnan(Cratio_ok)]
    #CratioC_ok = CratioC_ok[~np.isnan(CratioC_ok)]
    print 'snr_ymax_',  np.max(Fo)
    print 'snr_ymin_',  np.min(Fo)
    snr_ymax_ = np.max(Fo)
    snr_ymin_ = np.min(Fo)
  #  snr_plot.y_range.start = -0.2

   # exp_plot.yaxis.axis_label='Integration time for SNR = '+str(want_snr.value)+' [hours]' 

    snr_plot.title.text = 'Object Spectrum: '+template.value
  
    snr_plot.y_range.end = max(Fo)*1.1
    snr_plot.y_range.start = min(Fo)*0.9
    


    

       
######################################
# SET UP ALL THE WIDGETS AND CALLBACKS 
######################################

source = ColumnDataSource(data=dict(value=[]))
source.on_change('data', update_data)
exptime  = Slider(title="Integration Time (hours)", value=24., start=1., end=1000.0, step=1.0, callback_policy='mouseup')
exptime.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
magnitude  = Slider(title="V Band Magnitude of Object", value=10., start=-10, end=40., step=1.0, callback_policy='mouseup')
magnitude.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
diameter  = Slider(title="Mirror Diameter (meters)", value = 10.0, start=0.5, end=50., step=0.5, callback_policy='mouseup') 
diameter.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
resolution  = Slider(title="Telescope  Resolution", value = 100., start=10, end=10000, step=5, callback_policy='mouseup') 
resolution.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
temperature  = Slider(title="Telescope Temperature (K)", value = 270.0, start=90.0, end=400., step=10., callback_policy='mouseup') 
temperature.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
ntherm  = Slider(title="Number of thermal surfaces:", value = 1, start=1, end=30., step=1, callback_policy='mouseup') 
ntherm.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
darkcurrent  = Slider(title="Dark current (counts/s)", value = 1e-4, start=1e-5, end=1e-3, step=1e-5, callback_policy='mouseup') 
darkcurrent.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
readnoise  = Slider(title="Read noise (counts/pixel)", value = 0.1, start=0.01, end=1, step=0.05, callback_policy='mouseup') 
readnoise.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
dtmax  = Slider(title="Maximum single exposure time (hours)", value = 1, start=0.1, end=10., step=0.5, callback_policy='mouseup') 
dtmax.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
want_snr  = Slider(title="Desired signal-to-noise ratio? (only used for exposure time plot)", value = 10, start=0.5, end=100., step=0.5, callback_policy='mouseup') 
want_snr.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
gain  = Slider(title="Detector gain:", value = 1, start=1, end=1000., step=5, callback_policy='mouseup') 
gain.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
#ground based choice
ground_based = Select(title="Simulate ground-based observation?", value="No", options=["No",  "Yes"])

#select menu for planet
template = Select(title="Object Spectrum", value="Sun", options=["Sun" ])


oo = column(children=[exptime, diameter, resolution, temperature, ntherm, dtmax, gain, ground_based]) 
pp = column(children=[template, magnitude]) 
qq = column(children=[instruction0, text_input, instruction1, format_button_group, instruction2, link_box])
#ii = column(children=[inner, outer,  dtmax])
#ee = column(children=[want_snr])

observation_tab = Panel(child=oo, title='Observation')
planet_tab = Panel(child=pp, title='Target')
#instrument_tab = Panel(child=ii, title='Instrumentation')
download_tab = Panel(child=qq, title='Download')

for w in [text_input]: 
    w.on_change('value', change_filename)
format_button_group.on_click(i_clicked_a_button)

for ww in [template]: 
    ww.on_change('value', update_data)

for gg in [ground_based]: 
    gg.on_change('value', update_data)


inputs = Tabs(tabs=[ planet_tab, observation_tab,  download_tab ])

curdoc().add_root(row(inputs, ptabs)) 
