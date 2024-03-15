#!/usr/bin/env python
"""
Created on Tue Apr 25 10:19:25 2017

@author: gkanarek
"""

from __future__ import (print_function, division, absolute_import, with_statement,
                        nested_scopes, generators)
                        
#### Base Protocol class

class Protocol(object):
    """
    A class for defining new persistence protocols; the only one currently
    implemented is JSON.
    
    Each protocol needs to include four functions:
        1. An encode function
        2. A decode function
        3. A save function
        4. A load function
    
    Each protocol also needs to store units along with values in some way,
    and must at least handle attributes which are arrays or dicts.
    """
    
    def save(self, *args, **kwargs):
        raise NotImplementedError("This protocol is incomplete!")
    def load(self, *args, **kwargs):
        raise NotImplementedError("This protocol is incomplete!")
    def encode(self, *args, **kwargs):
        raise NotImplementedError("This protocol is incomplete!")
    def decode(self, *args, **kwargs):
        raise NotImplementedError("This protocol is incomplete!")

    
    def __init__(self, defaults_dict):
        self._defaults = defaults_dict