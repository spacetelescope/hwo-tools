import bokeh
from bokeh.models import ColumnDataSource, Div
from bokeh.plotting import figure, show
from bokeh.layouts import layout

import numpy as np

wave = np.arange(0.3,1.7,0.01)
flux = np.random.rand(len(wave))
sourcename = "Qualvo"

source = ColumnDataSource(data={"wave": wave, "flux": flux})

div = Div(text=r"<p>A sample spectrum received while scanning an <i>otherwise-unidentified sector</i> of space.</p>")
p = figure(width=750, height=300, title=f"Spectrum of {sourcename}", x_axis_label=r"$$\mu m$$", y_axis_label="Foobar", tools=("hover", "box_zoom", "wheel_zoom", "reset"), tooltips=[("@wave", "@flux")], toolbar_location="below")
p.line("wave", "flux", source=source, line_width=5, color="#921a35", legend_label="PHOTLAM")

design = layout([div], [p], sizing_mode="scale_width")

show(design)
