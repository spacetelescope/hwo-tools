# Import some standard python packages
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits, ascii 
from matplotlib import gridspec
from matplotlib import rc
import pdb
import sys
import os 
from astropy.table import Table, Column
mpl.rc('font', family='Times New Roman')
mpl.rcParams['font.size'] = 25.0

from bokeh.themes import Theme 
import yaml 
from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, HBox, VBoxForm, HoverTool, Paragraph, Range1d 
from bokeh.layouts import column, row, WidgetBox 
from bokeh.models.widgets import Slider, Panel, Tabs, Div, TextInput, RadioButtonGroup
from bokeh.io import hplot, vplot, curdoc
from bokeh.models.callbacks import CustomJS
from bokeh.embed import components, autoload_server 

import coronagraph as cg  # Import coronagraph model

################################
# PARAMETERS
################################

# Integration time (hours)
Dt = 20.0 # - SLIDER 

# Planet params
alpha = 90.     # phase angle at quadrature
Phi   = 1.      # phase function at quadrature (already included in SMART run)
Rp    = 1.0     # Earth radii - SLIDER 
r     = 1.0     # semi-major axis (AU) - SLIDER 

# Stellar params
Teff  = 5780.   # Sun-like Teff (K)
Rs    = 1.      # star radius in solar radii

# Planetary system params
d    = 10.     # distance to system (pc)  - SLIDER 
Nez  = 1.      # number of exo-zodis  - SLIDER 

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
    cg.count_rates(Ahr, lamhr, solhr, alpha, Phi, Rp, Teff, Rs, r, d, Nez)
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

planet = ColumnDataSource(data=dict(lam=lam, cratio=Cratio*1e9, spec=spec*1e9, downerr=(spec-sig)*1e9, uperr=(spec+sig)*1e9))

################################
# BOKEH PLOTTING
################################

snr_plot = Figure(plot_height=400, plot_width=750, 
              tools="crosshair,pan,reset,resize,save,box_zoom,wheel_zoom",
              x_range=[0.3, 2.7], y_range=[0, 1], toolbar_location='right')
snr_plot.x_range = Range1d(0.3, 2.7, bounds=(0.3, 2.7)) 
snr_plot.y_range = Range1d(0.0, 1.2, bounds=(0.3, 5.0)) 
snr_plot.background_fill_color = "beige"
snr_plot.background_fill_alpha = 0.5
snr_plot.yaxis.axis_label='F_p/F_s (x10^9)' 
snr_plot.xaxis.axis_label='Wavelength [micron]' 

snr_plot.line('lam','cratio',source=planet,line_width=2.0, color="green", alpha=0.7)
snr_plot.circle('lam', 'spec', source=planet, fill_color='red', line_color='black', size=8) 
snr_plot.segment('lam', 'downerr', 'lam', 'uperr', source=planet, line_width=1, line_color='grey', line_alpha=0.5) 

def change_filename(attrname, old, new): 
   format_button_group.active = None 


instruction0 = Div(text="""Choose a file rootname here 
                           (no special characters):""", width=300, height=15)
text_input = TextInput(value="filename", title=" ", width=100)
instruction1 = Div(text="""Then choose a format here:""", width=300, height=15)
format_button_group = RadioButtonGroup(labels=["txt", "fits"])
instruction2 = Div(text="""Your file will be linked here:""", width=300, height=15)
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



def update_data(attrname, old, new):

    print 'Updating model for exptime = ', exptime.value, ' for planet with R = ', radius.value, ' at distance ', distance.value, ' parsec '
    print '                   exozodi = ', exozodi.value
    
    # Run coronagraph with default LUVOIR telescope (aka no keyword arguments)
    lam, dlam, A, q, Cratio, cp, csp, cz, cez, cD, cR, cth, DtSNR = \
        cg.count_rates(Ahr, lamhr, solhr, alpha, Phi, radius.value, Teff, Rs, semimajor.value, distance.value, exozodi.value)
    # Calculate background photon count rates
    cb = (cz + cez + csp + cD + cR + cth)
    # Convert hours to seconds
    Dts = exptime.value * 3600.
    # Calculate signal-to-noise assuming background subtraction (the "2")
    SNR  = cp*Dts/np.sqrt((cp + 2*cb)*Dts)
    # Calculate 1-sigma errors
    sig= Cratio/SNR
    # Add gaussian noise to flux ratio
    spec = Cratio + np.random.randn(len(Cratio))*sig
    
    planet.data = dict(lam=lam, cratio=Cratio*1e9, spec=spec*1e9, downerr=(spec-sig)*1e9, uperr=(spec+sig)*1e9) 

    format_button_group.active = None 

######################################
# SET UP ALL THE WIDGETS AND CALLBACKS 
######################################


source = ColumnDataSource(data=dict(value=[]))
source.on_change('data', update_data)
exptime  = Slider(title="Integration Time (hours)", value=20., start=10., end=100.0, step=1.0, callback_policy='mouseup')
exptime.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
distance = Slider(title="Distance (parsec)", value=10., start=2., end=50.0, step=1.0, callback_policy='mouseup') 
distance.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
radius   = Slider(title="Planet Radius (R_Earth)", value=1.0, start=0.5, end=3., step=0.1, callback_policy='mouseup') 
radius.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
semimajor= Slider(title="Semi-major axis of orbit (AU)", value=1.0, start=0.2, end=2., step=0.1, callback_policy='mouseup') 
semimajor.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")
exozodi  = Slider(title="Number of Exozodi", value = 1.0, start=1.0, end=10., step=1., callback_policy='mouseup') 
exozodi.callback = CustomJS(args=dict(source=source), code="""
    source.data = { value: [cb_obj.value] }
""")




oo = column(children=[exptime]) 
pp = column(children=[distance, radius, semimajor, exozodi]) 
qq = column(children=[instruction0, text_input, instruction1, format_button_group, instruction2, link_box]) 

observation_tab = Panel(child=oo, title='Observation')
planet_tab = Panel(child=pp, title='Planet')
download_tab = Panel(child=qq, title='Download')

for w in [text_input]: 
    w.on_change('value', change_filename)
format_button_group.on_click(i_clicked_a_button) 


inputs = Tabs(tabs=[ planet_tab, observation_tab, download_tab ])
curdoc().add_root(row(children=[inputs, snr_plot])) 

#curdoc().theme = Theme(json=yaml.load("""
#attrs:
#    Figure:
#        background_fill_color: '#2F2F2F'
#        border_fill_color: '#2F2F2F'
#        outline_line_color: '#444444'
#    Axis:
#        axis_line_color: "white"
#        axis_label_text_color: "white"
#        major_label_text_color: "green"
#        major_tick_line_color: "white"
#        minor_tick_line_color: "white"
#        minor_tick_line_color: "white"
#    Grid:
#        grid_line_dash: [6, 4]
#        grid_line_alpha: .9
#    Title:
#        text_color: "green"
#""")) 


curdoc().add_root(source) 


