import numpy as np

class Output(object):

    def __init__(self, lam=None, dlam=None, A=None, q=None, Cratio=None,
                 cp=None, csp=None, cz=None, cez=None, cD=None, cR=None,
                 cth=None, DtSNR=None):

         self.lam = lam
         self.dlam = dlam
         self.A = A
         self.q = q
         self.Cratio = Cratio
         self.cp = cp
         self.csp = csp
         self.cz = cz
         self.cez = cez
         self.cD = cD
         self.cR = cR
         self.cth = cth
         self.DtSNR = DtSNR
