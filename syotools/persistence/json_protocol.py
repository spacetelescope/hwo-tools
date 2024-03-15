#!/usr/bin/env python
"""
Created on Tue Apr 25 11:43:22 2017

@author: gkanarek
"""

from __future__ import (print_function, division, absolute_import, 
                        with_statement, nested_scopes, generators)

#Protocols:
import json
from collections import OrderedDict

#Base class:
from .protocol import Protocol
from ..utils import JsonUnit, JsonSpectrum

####JSON protocol implementation
# 2017-11-20: converted to OrderedDict, to preserve attribute order.

class JSON(Protocol):
    
    name = 'json'
    
    def encode(self, entry):
        """
        Encode a single model attribute as a tuple or dict, for JSON storage.
        
        We'll use JsonUnit as the means of encoding/decoding astropy.unit.Quantity
        objects.
        """
        if isinstance(entry, (JsonUnit, JsonSpectrum)):
            return entry.encode_json()
        elif isinstance(entry, dict):
            return OrderedDict([(k, self.encode(v)) for k,v in entry.items()])
        else:
            return entry
    
    def decode(self, attr, entry):
        """
        Decode an individual entry in a JSON profile into the appropriate 
        Quantity(s). This method uses the 'generic' units formatting, to allow 
        for user-defined profiles that don't use the FITS standard.
        
        Special cases:
          - If the entry is not a dict or a tuple, it is returned without
            modification.
          - If the entry is a dict, recurse on each element.
          - If the entry is a list, convert it to a numpy array.
          - If the units can't be parsed, use the default unit for that
            attribute if one is available, otherwise use 
            astropy.units.UnrecognizedUnit.
        """
        if not isinstance(entry, (dict, tuple, OrderedDict)):
            return entry #This doesn't have units.
        if isinstance(entry, (dict, OrderedDict)):
            return OrderedDict([(k, self.decode(k, v)) for k,v in entry.items()])
        if "JsonUnit" in entry:
            return JsonUnit.decode_json(entry)
    
    def load(self, model_class, source):
        """
        Load a model profile from a json file (by default) or string. This is 
        a very basic persistence model, which can absolutely be made more 
        sophisticated.
        
        The only difference from standard JSON is that we want to parse each
        attribute through the protocol decoder.
        """
    
        try:
            with open(source) as f:
                profile_dict = json.load(f, object_pairs_hook=OrderedDict)
        except FileNotFoundError:
            print("File {} not found; using default profile.".format(source))
        
        return self.create_from_dict(profile_dict)
            
    def save(self, model, destination):
        """
        Store a model profile into a json file. To handle attributes with 
        units, we'll pass everything through the encoder.
        
        This is a very basic persistence model, which can absolutely be made 
        more sophisticated.
        """
        
        profile = self.encode_to_dict(model)

        with open(destination, 'w') as f:
            json.dump(f, profile, sort_keys=True, indent=4)
    
    def encode_to_dict(self, model):
        """
        Store a model profile as a dictionary, which can be saved in a file.
        """
        profile = {}

        for attr in model._tracked_attributes:
            val = getattr(model, attr)
            profile[attr] = self.encode(val)
        
        return profile
    
    def create_from_dict(self, model_class, profile_dict):
        """
        Create a model from a profile dictionary (presumably stored in a file).
        """
        profile = {}
        
        try:
            #Reconstruct the profile
            for attr, entry in profile_dict.items():
                profile[attr] = self.decode(attr, entry)
        except NotImplementedError:
            print("This protocol is not implemented; using default profile.")
            profile = OrderedDict()
        except json.JSONDecodeError:
            print("Unable to decode source; using default profile.")
            profile = OrderedDict()
        
        return model_class(**profile)
            
        
        