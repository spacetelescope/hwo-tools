#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 11:54:54 2017

@author: gkanarek

YAML I/O with OrderedDict to preserve element orders.

Note that this functionality is particularly important for the Exposure 
model, where we want to load _sed_id after _sed so that the default is 
loaded from pysynphot and not overwritten by the placeholder.
"""

from __future__ import (print_function, division, absolute_import, with_statement,
                        nested_scopes, generators)

import yaml
from collections import OrderedDict


#Loading and dumping functions to make use of OrderedDicts
def _construct_mapping(loader, node):
    loader.flatten_mapping(node)
    return OrderedDict(loader.construct_pairs(node))

def _dict_representer(dumper, data):
    return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
            data.items())


#yaml Loader and Dumper subclasses to use OrderedDicts
class OrderedLoader(yaml.Loader):
    pass

class OrderedDumper(yaml.Dumper):
    pass


#Register the loading & dumping functions  
OrderedLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _construct_mapping)

OrderedDumper.add_representer(OrderedDict, _dict_representer)


#Wrapper functions for yaml.load and yaml.dump to make use of OrderedLoader
#and OrderedDumper.

def ordered_load(stream):
    """
    Wrapper for yaml.load to make use of an OrderedDict. 
    """
    return yaml.load(stream, OrderedLoader)

def ordered_dump(data, stream=None, **kwds):
    """
    YAML serialization using an OrderedDict.
    """
    return yaml.dump(data, stream, OrderedDumper, **kwds)