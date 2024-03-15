#!/usr/bin/env python
"""
Created on Tue Apr 25 10:20:30 2017

@author: gkanarek
"""
from ..utils import DotDict
from .protocol import Protocol

from .json_protocol import JSON

protocols = DotDict({'json':JSON})
