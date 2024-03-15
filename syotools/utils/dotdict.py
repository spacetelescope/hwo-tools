#!/usr/bin/env python
"""
Created on Tue Apr 25 10:24:25 2017

@author: gkanarek

Utility class for dot access to a dict. Adapted from:
http://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary
"""

class DotDict(dict):
    """
    A subclass of dict to allow object dot access.
    """

    def __getattr__(self, attr):
        if attr not in self:
            nm = self.__class__.__name__
            raise AttributeError("{} object has no attribute '{}'".format(nm, attr))
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(DotDict, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(DotDict, self).__delitem__(key)
        del self.__dict__[key]