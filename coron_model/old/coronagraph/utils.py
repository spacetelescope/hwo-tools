import imp, sys
from types import ModuleType, FunctionType, StringType
import os

__all__ = ['Input']

inpath = "inputs/"
relpath = os.path.join(os.path.dirname(__file__), inpath)

class Input(object):
    """
    Reads default and user input files and creates a class where the input file
    variables are attributes. See inputs/ to customize.

    """
    def __init__(self, input_type = ''):

        if (input_type == 'telescope') or (input_type == 'planet') or (input_type == 'star'):
            pass
        else:
            print "Error: unrecognized input_type. Please use 'telescope', 'planet', or 'star'."
            return

        try:
            del sys.modules['input_usr']
            del sys.modules['input']
        except KeyError:
            pass

        default_input_file = os.path.join(relpath,'input_default_'+input_type+'.py')
        user_input_file = os.path.join(relpath,'input_user_'+input_type+'.py')

        self._input = imp.load_source("input", default_input_file)            # Load default inputs into self._input

        self._input_usr = imp.load_source("input_usr", user_input_file)       # Load user inputs into self._input_usr

        self._input.__dict__.update(self._input_usr.__dict__)                 # Update self._input with user values

        inp_dict = self._input.__dict__

        for key, value in inp_dict.items():
            if key.startswith('__') or isinstance(value, ModuleType) or isinstance(value, FunctionType):
                inp_dict.pop(key, None)

        self.__dict__.update(inp_dict)                                        # Make all parameters accessible as self.param

        del self._input
        del self._input_usr

class Loadin(object):
    """
    Reads default and user input files and creates a class where the input file
    variables are attributes. See inputs/ to customize.

    """
    def __init__(self, path):

        if not path.endswith('.py'):
            print("Incompatible file.")
            return

        try:
            del sys.modules['input']
        except KeyError:
            pass

        user_input_file = path

        self._input = imp.load_source("input", user_input_file)            # Load inputs into self._input

        inp_dict = self._input.__dict__

        for key, value in inp_dict.items():
            if key.startswith('__') or isinstance(value, ModuleType) or isinstance(value, FunctionType):
                inp_dict.pop(key, None)

        self.__dict__.update(inp_dict)                                        # Make all parameters accessible as self.param

        del self._input
