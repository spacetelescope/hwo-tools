#!/usr/bin/env python2
"""
Created on Mon Feb 20 14:05:03 2017

@author: gkanarek
"""

from bokeh.layouts import column, row, widgetbox, gridplot, layout
from bokeh.plotting import figure as Figure
from bokeh.models.tools import HoverTool
from bokeh.models.ranges import Range1d
from bokeh.models.sources import ColumnDataSource
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import (Slider, Tabs, Div, Panel, Select, TextInput,
                                  Button, Toggle, RadioButtonGroup)
from bokeh.io import curdoc

mappings = {'CustomJS': CustomJS,
            'Range1d': Range1d,
            'ColumnDataSource': ColumnDataSource,
            'HoverTool': HoverTool,
            'Slider': Slider,
            'Panel': Panel,
            'Div': Div,
            'Select': Select,
            'Tabs': Tabs,
            'widgetbox': widgetbox,
            'TextInput': TextInput,
            'Button': Button,
            'gridplot': gridplot,
            'Toggle': Toggle,
            'RadioButtonGroup': RadioButtonGroup,
            'layout': layout
            }

sequences = {'column': column,
             'row': row,
            }

scalars = {}

def mapping_factory(tool, element_type):
    def mapping_constructor(loader, node):
        fmt = tool.formats.get(element_type, {})
        value = loader.construct_mapping(node, deep=True)
        ref = value.pop("ref", "")
        callback = value.pop("on_change", [])
        onclick = value.pop("on_click", None)
        fmt.update(value)
        if element_type == "Slider":
            fmt["start"], fmt["end"], fmt["step"] = fmt.pop("range", [0, 1, 0.1])
        arg = fmt.pop("arg", None)
        if arg is not None:
            obj = mappings[element_type](arg, **fmt)
        else:
            obj = mappings[element_type](**fmt)
        if ref:
            tool.refs[ref] = obj
        if callback:
            obj.on_change(*callback)
        if onclick:
            obj.on_click(onclick)
        yield obj
    
    mapping_constructor.__name__ = element_type.lower() + '_' + mapping_constructor.__name__
    return mapping_constructor

def sequence_factory(tool, element_type):
    def sequence_constructor(loader, node):
        fmt = tool.formats.get(element_type, {})
        value = loader.construct_sequence(node, deep=True)
        #ref = value.pop("ref", "") #can't have these in a sequence
        #callback = value.pop("on_change", []) #can't have these in a sequence
        obj = sequences[element_type](*value, **fmt)
        #if ref:
        #    tool.refs[ref] = obj
        #if callback:
        #    obj.on_change(*callback)
        yield obj
        
    sequence_constructor.__name__ = element_type.lower() + '_' + sequence_constructor.__name__
    return sequence_constructor

def scalar_factory(tool, element_type):
    def scalar_constructor(loader, node):
        fmt = tool.formats.get(element_type, {})
        value = loader.construct_scalar(node)
        ref = value.pop("ref", "")
        callback = value.pop("on_change", [])
        obj = scalars[element_type](value, **fmt)
        if ref:
            tool.refs[ref] = obj
        if callback:
            obj.on_change(*callback)
        yield obj
        
    scalar_constructor.__name__ = element_type.lower() + '_' + scalar_constructor.__name__
    return scalar_constructor

#These constructors need more specialized treatment

def document_constructor(tool, loader, node):
    layout = loader.construct_sequence(node, deep=True)
    for element in layout:
        curdoc().add_root(element)
    tool.document = curdoc()
    yield tool.document

def figure_constructor(tool, loader, node):
    
    fig = loader.construct_mapping(node, deep=True)
    fmt = tool.formats.get('Figure', {})
    
    elements = fig.pop('elements', [])
    cmds = []
    ref = fig.pop("ref", "")
    callback = fig.pop("on_change", [])
    axis = tool.formats.get("Axis", {})
    
    for key in fig:
        val = fig[key]
        if key in ['text', 'add_tools']:
           cmds.append((key, val))
        else:
            fmt[key] = val
    
    figure = Figure(**fmt)
    
    for key, cmd in cmds:
        if key == 'add_tools':
            figure.add_tools(*cmd)
        elif key == 'text':
            figure.text(*cmd.pop('loc'), **cmd)
    
    for element in elements:
        key = element.pop('kind')
        if key == 'line':
            line_fmt = tool.formats.get('Line', {})
            line_fmt.update(element)
            figure.line('x', 'y', **line_fmt)
        elif key == 'circle':
            circle_fmt = tool.formats.get('Circle', {})
            circle_fmt.update(element)
            figure.circle('x', 'y', **circle_fmt)
        elif key == 'image':
            image_fmt = tool.formats.get('Image', {})
            image_fmt.update(element)
            arg = image_fmt.pop("image", None)
            figure.image_rgba(arg, **image_fmt)
            
    for attr, val in axis.items():
        #change axis attributes, hopefully
        setattr(figure.axis, attr, val)
    
    if ref:
        tool.refs[ref] = figure
    if callback:
        figure.on_change(*callback)

    yield figure
        