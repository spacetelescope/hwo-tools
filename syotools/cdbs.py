"""
Created on Wed Mar 29 12:02:37 2017

@author: gkanarek
"""

import os

if 'PYSYN_CDBS' not in os.environ:
    os.environ['PYSYN_CDBS'] = os.path.expanduser("~/cdbs") 