from __future__ import (print_function, division, absolute_import, with_statement,
                        nested_scopes, generators)

"""try:
    import .hv_factory
    _HV = True
except ImportError:
    from warnings import warn
    warn("holoviews not installed - only standard Bokeh currently available",
         ImportWarning)
    import .std_factory
    _HV = False

factory = hv_factory if _HV else std_factory"""

from .base import SYOTool
