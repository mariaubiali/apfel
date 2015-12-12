************************************************************************
*
*     CleanUp.f:
*
*     It unsets all the evolution parameters so that they need to be 
*     reset.
*
************************************************************************
      subroutine CleanUp
*
      implicit none
*
      include "../commons/InAPFEL.h"
      include "../commons/InAPFELDIS.h"
      include "../commons/Welcome.h"
      include "../commons/scales.h"
      include "../commons/Evs.h"
      include "../commons/Nf_FF.h"
      include "../commons/ipt.h"
      include "../commons/Th.h"
      include "../commons/alpha_ref_QCD.h"
      include "../commons/alpha_ref_QED.h"
      include "../commons/kren.h"
      include "../commons/mass_scheme.h"
      include "../commons/m2th.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/grid.h"
      include "../commons/pdfset.h"
      include "../commons/Replica.h"
      include "../commons/lock.h"
      include "../commons/EvolOp.h"
      include "../commons/lambda_ref_QCD.h"
      include "../commons/AlphaEvolution.h"
      include "../commons/PDFEvolution.h"
      include "../commons/TimeLike.h"
      include "../commons/Polarized.h"
      include "../commons/Smallx.h"
      include "../commons/FastEvol.h"
      include "../commons/MassScheme.h"
      include "../commons/ProcessDIS.h"
      include "../commons/PolarizationDIS.h"
      include "../commons/ProjectileDIS.h"
      include "../commons/TargetDIS.h"
      include "../commons/ZedMass.h"
      include "../commons/WMass.h"
      include "../commons/ProtonMass.h"
      include "../commons/Sin2ThetaW.h"
      include "../commons/GFermi.h"
      include "../commons/CKM.h"
      include "../commons/MassRunning.h"
      include "../commons/TMC.h"
      include "../commons/DampingFONLL.h"
      include "../commons/SelectedCharge.h"
      include "../commons/TauMass.h"
      include "../commons/LeptEvol.h"
      include "../commons/LHAgrid.h"
      include "../commons/krenQ.h"
      include "../commons/kfacQ.h"
      include "../commons/PropagatorCorrection.h"
      include "../commons/EWCouplings.h"
      include "../commons/DynScVar.h"
      include "../commons/IntrinsicCharm.h"
*
*     Set all the initialization flags to "xxxx" so that
*     they will be initialized again by InitializeAPFEL.
*
      InAPFEL           = "xxxx"
*
      InWelcome         = "xxxx"
      InScales          = "xxxx"
      InPt              = "xxxx"
      InEvs             = "xxxx"   
      InTheory          = "xxxx"
      InAlpQCD          = "xxxx"
      InAlpQED          = "xxxx"
      InKren            = "xxxx"  
      InMasses          = "xxxx"
      InMassRef         = "xxxx"
      InThrRatios       = "xxxx"
      InMTau            = "xxxx"
      InMassRunning     = "xxxx"
      InMFP             = "xxxx"
      InMFA             = "xxxx"
      InPDFs            = "xxxx"
      InRep             = "xxxx"
      InEvolOp          = "xxxx"
      InLeptEvol        = "xxxx"
      InLock            = "xxxx"
      InGrid            = "xxxx"
      InTimeLike        = "xxxx"
      InPolarized       = "xxxx"
      InSmallx          = "xxxx"
      InAlphaEvol       = "xxxx"
      InLambdaQCD       = "xxxx"
      InPDFEvol         = "xxxx"
      InFastEvol        = "xxxx"
      InLHgrid          = "xxxx"
*
      InAPFELDIS        = "xxxx"
*
      InMassScheme      = "xxxx"
      InProcessDIS      = "xxxx"
      InPolarizationDIS = "xxxx"
      InProjectileDIS   = "xxxx"
      InTargetDIS       = "xxxx"
      InTMC             = "xxxx"
      InDampingFONLL    = "xxxx"
      InSelectedCharge  = "xxxx"
      InKrenQ           = "xxxx"
      InKfacQ           = "xxxx"
      InDynScVar        = "xxxx"
      InIntrinsicCharm  = "xxxx"
*
      InMZ              = "xxxx"
      InMW              = "xxxx"
      InMProton         = "xxxx"
      InSin2ThetaW      = "xxxx"
      InGFermi          = "xxxx"
      InCKM             = "xxxx"
      InDeltaR          = "xxxx"
      InEWCouplings     = "xxxx"
*
      return
      end
