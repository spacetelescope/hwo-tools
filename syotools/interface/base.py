#!/usr/bin/env python2
"""
Created on Wed Oct 19 12:43:00 2016

@author: gkanarek, tumlinson
"""
from __future__ import (print_function, division, absolute_import, 
                        with_statement, nested_scopes, generators)

import yaml
import os
import json
from random import sample
from string import digits, ascii_letters

#Import the constructor function factories
from syotools.interface import factory


class SYOParserError(Exception):
    """
    A custom error for problems with parsing the interface files.
    """


class SYOTool(object):
    """
    This is a framework for quickly generating Bokeh tools for Science Yield
    Optimization. Much of the Bokeh machinery will be kept under the hood in
    this class. To create a new tool, simply subclass SYOTool, set some
    parameters, connect with some astronomical models, and output embeddable 
    HTML.
    
    One major update is to use .yaml files as a templating method (a la the
    kv language in Kivy) to allow for quick and easy human-readable Bokeh
    interface construction.
    """
    
    #Here we assign the constructor methods
    mapping_factory = factory.mapping_factory
    sequence_factory = factory.sequence_factory
    scalar_factory = factory.scalar_factory
    figure_constructor = factory.figure_constructor
    document_constructor = factory.document_constructor
    
    hvfigure_constructor = factory.hvfigure_constructor
    renderer_constructor = factory.renderer_constructor
    
    #Attributes required for saving calculations
    tool_prefix = None #must be set by the subclass
    save_ext = '.json'
    save_dir = "saves"
    save_models = [] #must be set by the subclass
    save_params = {} #must be set by the subclass
    tool_defaults = {}
    load_mismatch = False #do the loaded models and parameters match?
    
    #storing the existing save JSON dictionary, in case of cross-tool calculations
    _save_json = {'models':{}}
    
        
    def self_constructor(self, loader, tag_suffix, node):
        """
        A multi_constructor for `!self` tag in the interface file.
        """
        yield eval("self"+tag_suffix, globals(), locals())
        
    def tool_preinit(self):
        """
        This should be implemented by the tool subclass, to do any pre-
        initialization steps that the tool requires.
        
        If this is not required, subclass should set `tool_preinit = None`
        in the class definition.
        """
        
        raise NotImplementedError
    
    def tool_postinit(self):
        """
        This should be implemented by the tool subclass, to do any pre-
        initialization steps that the tool requires.
        
        If this is not required, subclass should set `tool_postinit = None`
        in the class definition.
        """
        
        raise NotImplementedError
    
    def __init__(self):
        """
        First step is to set all the tool default parameters.
        
        Then, we need to initialize two dictionary attributes:
            
            self.formats, to include any formatting keywords which are set 
                in a separate string
                
            self.refs, which contains all the specific Bokeh objects which
                are created by the YAML parser.
        
        Then, we need to register all the constructors that we need.
        """
        
        #Parse tool parameter defaults.
        for param, val in self.tool_defaults.items():
            setattr(self, param, val)
        
        #Allow for pre-init stuff from the tool subclass.
        if self.tool_preinit is not None:
            self.tool_preinit()
        
        
        #Handle YAML construction
        self.formats = {}
        self.refs = {}
        
        self.document = None
        
        #Register constructors
        for m in factory.mappings:
            yaml.add_constructor(u"!"+m+":", self.mapping_factory(m))
        for s in factory.sequences:
            yaml.add_constructor(u"!"+s+":", self.sequence_factory(s))
        for s in factory.scalars:
            print(s)
            yaml.add_constructor(u"!"+s+":", self.scalar_factory(s))
        
        yaml.add_constructor(u"!Figure:", self.figure_constructor)
        yaml.add_constructor(u"!Document:", self.document_constructor)
        yaml.add_constructor(u"!renderer:", self.renderer_constructor)
        yaml.add_constructor(u"!HVFigure:", self.hvfigure_constructor)
        
        yaml.add_multi_constructor(u"!self", self.self_constructor)
        
        self.include_formatting()
        self.parse_interface()
        
        #Allow for post-init stuff from the tool subclass.
        if self.tool_postinit is not None:
            self.tool_postinit()
    
    def include_formatting(self):
        """
        This should simply be a dictionary of formatting keywords.
        """
        if not self.format_string:
            return
        
        self.formats = yaml.load(self.format_string)
    
    def parse_interface(self):
        """
        This is the workhorse YAML parser, which creates the interface based
        on the layout file.
        
        `interface_file` is the path to the interface .yaml file to be parsed.
        """
        
        if not self.interface_file:
            raise NotImplementedError("Interface file required.")
        
        #Read the interface file into a string
        filepath = os.path.abspath(os.path.expanduser(self.interface_file))
        if not os.path.exists(filepath):
            raise SYOParserError("Interface file path does not exist.")
        with open(filepath) as f:
            interface = f.read()
        
        #First, let's make sure that there's a Document in here
        if '!Document' not in interface:
            raise SYOParserError("Interface file must contain a Document tag")
        
        #Now, since we've registered all the constructors, we can parse the 
        #entire string with yaml. We don't need to assign the result to a
        #variable, since the constructors store everything in self.refs
        #(and self.document, for the document)
        
        self.full_stream = list(yaml.load_all(interface))
    
    def validate_filename(self, savefile="", overwrite=False):
        if not self.tool_prefix:
            raise NotImplementedError("Tool prefix identifier required to save calculations")
        
        if not hasattr(self, 'user_prefix') or not self.user_prefix:
            raise NotImplementedError("User prefix identifier required to save calculations")
        
        if not savefile:
            postfix = ''.join(sample(digits + ascii_letters, 8))
            savefile = "{}.{}".format(self.user_prefix, postfix)
            testfile = savefile + self.save_ext
            
        while not overwrite and os.path.exists(os.path.join(self.save_dir, 
                                                            testfile)):
            postfix = sample(digits + ascii_letters, 8)
            savefile = "{}.{}".format(self.user_prefix, postfix)
            testfile = savefile + self.save_ext
        
        return savefile
    
    def save_file(self, savefile="", overwrite=False):
        """
        Save the relevant models and parameters to a file.
        
        As of 2017-10-30, this is being reworked to implement the following
        save/load scheme:
            
        - Top-level dictionary split into two sub-dictionaries: 
            - "tools", including a dictionary for the parameters associated 
              with each individual SYOTool subclass which has used the same 
              save file;
            - "models", including a dictionary for the parameters associated 
              with the individual model instances involved in a calculation.
        
        When saving from a particular tool instance, the model parameters will
        be updated (overwriting the previous parameters, if any), and the tool
        parameters FOR THAT TOOL ONLY will be updated.
        
        An individual calculation is identified by its exposure id (stored in
        the Exposure.exp_id parameter), which is also used to construct the
        filename. Whenever a new calculation file is saved, a new exp_id is
        generated for the Exposure model.
        """
            
        if (not self.save_params and not self.save_models):
            raise NotImplementedError("No models or parameters to save.")
            
        #Collect parameters and models we want to save
        par_key = '{}_params'.format(self.tool_prefix)
        self._save_json[par_key] = {param: getattr(self, param) \
                                       for param in self.save_params}
        
        #Only update models that are tracked by this tool
        for model in self.save_models:
            self._save_json['models'][model] = getattr(self, model).encode()
        
        #Make sure we have a workable filename
        filename = self.validate_filename(savefile, overwrite=overwrite)
        outfile = filename + self.save_ext
        
        try:
            #Dump to the indicated file
            with open(os.path.join(self.save_dir, outfile), 'w') as f:
                json.dump(self._save_json, f)
        except Exception as e:
            print(e)
            return ""
        
        self.load_mismatch = False #by definition, it now matches
        return filename
    
    def _check_param(self, param_name, model_name, model_attr):
        """
        Verify that the value of the indicated parameter matches the value
        stored in the mapped model attribute.
        """
        
        par = getattr(self, param_name)
        model = getattr(self, model_name)
        
        return par == getattr(model, model_attr)
    
    def load_file(self, savefile):
        """
        Load the stored calculation's models and parameters. Returns 0 on a 
        successful load; 1 means the savefile doesn't exist, 2 means there was
        an error retrieving the save state; 3 means an error in loading param
        and model values.
        
        As of 2017-10-30, this is being reworked to implement the following
        save/load scheme:
            
        - Top-level dictionary split into two sub-dictionaries: 
            - "tools", including a dictionary for the parameters associated 
              with each individual SYOTool subclass which has used the same 
              save file;
            - "models", including a dictionary for the parameters associated 
              with the individual model instances involved in a calculation.
        
        When loading into a particular tool instance, the model parameters will
        be loaded, and will be used to set the parameters of the tool, if 
        possible. Three notes:
            - Any tool parameters which are not set by model parameters will
              revert to the values in the save dict (as long as they don't
              cause some sort of conflict).
            - If there is no save dictionary for the current tool included in
              the save file, a warning will be provided to the user, and any 
              model-independent parameters will be set to the default values 
              for that tool.
            - If there is a save dictionary for the current tool included in
              the save file, and if the model-dependent tool parameters do not
              match the values in their respective model dictionaries, then
              a warning will be provided to the user, and the model values will
              take priority.
        """
        
        if (not self.save_params and not self.save_models):
            raise NotImplementedError("No models or parameters to load.")
        
        loadfile = os.path.join(self.save_dir, savefile + self.save_ext)
        if not os.path.exists(loadfile):
            return 1
        
        try:
            #Load the savestate
            with open(loadfile) as f:
                state = json.load(f)
        except Exception as e:
            print(e)
            return 2
        try:
            #Restore tracked parameters and models from the savestate (leave alone
            #anything which is not included in the savestate). Also check for
            #consistency between any save_params mapped to particular model 
            #attributes
            for model in self.save_models:
                if model in state['models']:
                    model_state = getattr(self, model)
                    model_state.decode(state['models'][model])
                    setattr(self, model, model_state)
            
            par_key = '{}_params'.format(self.tool_prefix)
            if par_key not in state:
                return 0
            for param, model_map in self.save_params.items():
                if param in state[par_key]:
                    setattr(self, param, state[par_key][param])
                    if not self.load_mismatch and model_map is not None:
                        self.load_mismatch = self._check_param(param, *model_map)
        except Exception as e:
            print(e)
            return 3
        return 0