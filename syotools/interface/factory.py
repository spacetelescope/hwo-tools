#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 17:48:52 2018

@author: gkanarek
"""

from __future__ import absolute_import


import holoviews as hv
from holoviews.operation.datashader import datashade
import parambokeh
from bokeh.io import curdoc

hv.extension('bokeh')
_renderer = hv.renderer('bokeh').instance(mode='server')

from ._std_factory import (mappings as std_mappings, 
                           sequences as std_sequences, 
                           scalars as std_scalars, 
                           mapping_factory as std_mapping_factory, 
                           scalar_factory as std_scalar_factory,
                           sequence_factory,
                           figure_constructor,
                           document_constructor)

mappings = std_mappings.copy()
sequences = std_sequences.copy()
scalars = std_scalars.copy()

_hv_mappings = {'HVPoints': hv.Points,
                'HVStream': hv.streams.Stream.define,
                'HVDynamic': hv.util.Dynamic,
                'PWidget': parambokeh.Widgets,
                'Datashader': datashade}

mappings.update(_hv_mappings)

_hv_scalars = {'HVRangeXY': hv.streams.RangeXY}

scalars.update(_hv_scalars)

def mapping_factory(tool, element_type):
    if element_type not in _hv_mappings:
        return std_mapping_factory(tool, element_type)
    def mapping_constructor(loader, node):
        fmt = tool.formats.get(element_type, {})
        value = loader.construct_mapping(node, deep=True)
        ref = value.pop("ref", "")
        callback = value.pop("on_change", [])
        onclick = value.pop("on_click", None)
        fmt.update(value)
        arg = fmt.pop("arg", None)
        if element_type == "HVStream":
            obj = mappings[element_type](arg, **fmt)()
        elif arg is not None:
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

def scalar_factory(tool, element_type):
    if element_type not in _hv_scalars:
        return std_scalar_factory(tool, element_type)
    def scalar_constructor(loader, node):
        fmt = tool.formats.get(element_type, {})
        value = loader.construct_scalar(node)
        #ref = value.pop("ref", "")
        #callback = value.pop("on_change", [])
        if not value:
            yield scalars[element_type]
        else:
            obj = scalars[element_type](value, **fmt)
            #if ref:
            #    tool.refs[ref] = obj
            #if callback:
            #    obj.on_change(*callback)
            yield obj
        
    scalar_constructor.__name__ = element_type.lower() + '_' + scalar_constructor.__name__
    return scalar_constructor

def renderer_constructor(tool, loader, node):
    opts = loader.construct_mapping(node, deep=True)
    ref = opts.pop("ref", "")
    opts.update(tool.formats.get("renderer", {}))
    for o in opts["opts"]:
        hv.opts(o)
    if ref:
        tool.refs[ref] = _renderer
    yield _renderer

def hvfigure_constructor(tool, loader, node):
    #Grab all of the formatting
    fig = loader.construct_mapping(node, deep=True)
    fmt = tool.formats.get('HVFigure', {})
    
    elements = fig.pop('elements', [])
    cmds = []
    ref = fig.pop("ref", "")
    callback = fig.pop("on_change", [])
    axis = tool.formats.get("Axis", {})
    xaxis = fig.pop("xaxis", {})
    yaxis = fig.pop("yaxis", {})

    data = fig.pop("source", None)
    
    for key in fig:
        val = fig[key]
        if key in ['text', 'add_tools']:
           cmds.append((key, val))
        else:
            fmt[key] = val
    
    #Create plot object and the figure itself
    plot_obj = _renderer.get_plot(data, doc=curdoc())
    figure = plot_obj.state
    
    #Since we can't just use the keyword dictionary with get_plot, we have to 
    #use setattr for all the formatting
    for key, val in fmt.items():
        setattr(figure, key, val)
    
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
        elif key == 'quad':
            quad_fmt = tool.formats.get('Quad', {})
            quad_fmt.update(element)
            figure.quad(**quad_fmt)
        
    for attr, val in xaxis.items() + axis.items():
        #change x-axis attributes, hopefully
        setattr(figure.xaxis, attr, val)
    
    for attr, val in yaxis.items() + axis.items():
        #change y-axis attributes, hopefully
        setattr(figure.yaxis, attr, val)
    
    if ref:
        tool.refs[ref] = figure
    if callback:
        figure.on_change(*callback)

    yield figure
