# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.11
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_apfel', [dirname(__file__)])
        except ImportError:
            import _apfel
            return _apfel
        if fp is not None:
            try:
                _mod = imp.load_module('_apfel', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _apfel = swig_import_helper()
    del swig_import_helper
else:
    import _apfel
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _apfel.delete_SwigPyIterator
    __del__ = lambda self : None;
    def value(self): return _apfel.SwigPyIterator_value(self)
    def incr(self, n=1): return _apfel.SwigPyIterator_incr(self, n)
    def decr(self, n=1): return _apfel.SwigPyIterator_decr(self, n)
    def distance(self, *args): return _apfel.SwigPyIterator_distance(self, *args)
    def equal(self, *args): return _apfel.SwigPyIterator_equal(self, *args)
    def copy(self): return _apfel.SwigPyIterator_copy(self)
    def next(self): return _apfel.SwigPyIterator_next(self)
    def __next__(self): return _apfel.SwigPyIterator___next__(self)
    def previous(self): return _apfel.SwigPyIterator_previous(self)
    def advance(self, *args): return _apfel.SwigPyIterator_advance(self, *args)
    def __eq__(self, *args): return _apfel.SwigPyIterator___eq__(self, *args)
    def __ne__(self, *args): return _apfel.SwigPyIterator___ne__(self, *args)
    def __iadd__(self, *args): return _apfel.SwigPyIterator___iadd__(self, *args)
    def __isub__(self, *args): return _apfel.SwigPyIterator___isub__(self, *args)
    def __add__(self, *args): return _apfel.SwigPyIterator___add__(self, *args)
    def __sub__(self, *args): return _apfel.SwigPyIterator___sub__(self, *args)
    def __iter__(self): return self
SwigPyIterator_swigregister = _apfel.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)


def InitializeAPFEL():
  return _apfel.InitializeAPFEL()
InitializeAPFEL = _apfel.InitializeAPFEL

def EvolveAPFEL(*args):
  return _apfel.EvolveAPFEL(*args)
EvolveAPFEL = _apfel.EvolveAPFEL

def xPDF(*args):
  return _apfel.xPDF(*args)
xPDF = _apfel.xPDF

def dxPDF(*args):
  return _apfel.dxPDF(*args)
dxPDF = _apfel.dxPDF

def xPDFj(*args):
  return _apfel.xPDFj(*args)
xPDFj = _apfel.xPDFj

def xgamma(*args):
  return _apfel.xgamma(*args)
xgamma = _apfel.xgamma

def xgammaj(*args):
  return _apfel.xgammaj(*args)
xgammaj = _apfel.xgammaj

def dxgamma(*args):
  return _apfel.dxgamma(*args)
dxgamma = _apfel.dxgamma

def ExternalEvolutionOperator(*args):
  return _apfel.ExternalEvolutionOperator(*args)
ExternalEvolutionOperator = _apfel.ExternalEvolutionOperator

def LHAPDFgrid(*args):
  return _apfel.LHAPDFgrid(*args)
LHAPDFgrid = _apfel.LHAPDFgrid

def LHAPDFgridDerivative(*args):
  return _apfel.LHAPDFgridDerivative(*args)
LHAPDFgridDerivative = _apfel.LHAPDFgridDerivative

def AlphaQCD(*args):
  return _apfel.AlphaQCD(*args)
AlphaQCD = _apfel.AlphaQCD

def AlphaQED(*args):
  return _apfel.AlphaQED(*args)
AlphaQED = _apfel.AlphaQED

def HeavyQuarkMass(*args):
  return _apfel.HeavyQuarkMass(*args)
HeavyQuarkMass = _apfel.HeavyQuarkMass

def NPDF(*args):
  return _apfel.NPDF(*args)
NPDF = _apfel.NPDF

def Ngamma(*args):
  return _apfel.Ngamma(*args)
Ngamma = _apfel.Ngamma

def LUMI(*args):
  return _apfel.LUMI(*args)
LUMI = _apfel.LUMI

def GetVersion():
  return _apfel.GetVersion()
GetVersion = _apfel.GetVersion

def CleanUp():
  return _apfel.CleanUp()
CleanUp = _apfel.CleanUp

def EnableWelcomeMessage(*args):
  return _apfel.EnableWelcomeMessage(*args)
EnableWelcomeMessage = _apfel.EnableWelcomeMessage

def EnableEvolutionOperator(*args):
  return _apfel.EnableEvolutionOperator(*args)
EnableEvolutionOperator = _apfel.EnableEvolutionOperator

def LockGrids(*args):
  return _apfel.LockGrids(*args)
LockGrids = _apfel.LockGrids

def SetTimeLikeEvolution(*args):
  return _apfel.SetTimeLikeEvolution(*args)
SetTimeLikeEvolution = _apfel.SetTimeLikeEvolution

def SetFastEvolution(*args):
  return _apfel.SetFastEvolution(*args)
SetFastEvolution = _apfel.SetFastEvolution

def SetSmallxResummation(*args):
  return _apfel.SetSmallxResummation(*args)
SetSmallxResummation = _apfel.SetSmallxResummation

def SetAlphaQCDRef(*args):
  return _apfel.SetAlphaQCDRef(*args)
SetAlphaQCDRef = _apfel.SetAlphaQCDRef

def SetAlphaQEDRef(*args):
  return _apfel.SetAlphaQEDRef(*args)
SetAlphaQEDRef = _apfel.SetAlphaQEDRef

def SetAlphaEvolution(*args):
  return _apfel.SetAlphaEvolution(*args)
SetAlphaEvolution = _apfel.SetAlphaEvolution

def SetLambdaQCDRef(*args):
  return _apfel.SetLambdaQCDRef(*args)
SetLambdaQCDRef = _apfel.SetLambdaQCDRef

def SetPDFEvolution(*args):
  return _apfel.SetPDFEvolution(*args)
SetPDFEvolution = _apfel.SetPDFEvolution

def SetQLimits(*args):
  return _apfel.SetQLimits(*args)
SetQLimits = _apfel.SetQLimits

def SetFFNS(*args):
  return _apfel.SetFFNS(*args)
SetFFNS = _apfel.SetFFNS

def SetGridParameters(*args):
  return _apfel.SetGridParameters(*args)
SetGridParameters = _apfel.SetGridParameters

def SetExternalGrid(*args):
  return _apfel.SetExternalGrid(*args)
SetExternalGrid = _apfel.SetExternalGrid

def SetMaxFlavourAlpha(*args):
  return _apfel.SetMaxFlavourAlpha(*args)
SetMaxFlavourAlpha = _apfel.SetMaxFlavourAlpha

def SetMaxFlavourPDFs(*args):
  return _apfel.SetMaxFlavourPDFs(*args)
SetMaxFlavourPDFs = _apfel.SetMaxFlavourPDFs

def SetMSbarMasses(*args):
  return _apfel.SetMSbarMasses(*args)
SetMSbarMasses = _apfel.SetMSbarMasses

def SetNumberOfGrids(*args):
  return _apfel.SetNumberOfGrids(*args)
SetNumberOfGrids = _apfel.SetNumberOfGrids

def SetPDFSet(*args):
  return _apfel.SetPDFSet(*args)
SetPDFSet = _apfel.SetPDFSet

def SetPerturbativeOrder(*args):
  return _apfel.SetPerturbativeOrder(*args)
SetPerturbativeOrder = _apfel.SetPerturbativeOrder

def SetPoleMasses(*args):
  return _apfel.SetPoleMasses(*args)
SetPoleMasses = _apfel.SetPoleMasses

def SetRenFacRatio(*args):
  return _apfel.SetRenFacRatio(*args)
SetRenFacRatio = _apfel.SetRenFacRatio

def SetReplica(*args):
  return _apfel.SetReplica(*args)
SetReplica = _apfel.SetReplica

def SetTheory(*args):
  return _apfel.SetTheory(*args)
SetTheory = _apfel.SetTheory

def SetVFNS():
  return _apfel.SetVFNS()
SetVFNS = _apfel.SetVFNS

def DIS_xsec(*args):
  return _apfel.DIS_xsec(*args)
DIS_xsec = _apfel.DIS_xsec
# This file is compatible with both classic and new-style classes.


