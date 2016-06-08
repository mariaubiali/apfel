************************************************************************
*
*     initHELL.f:
*
*     This routine initializes the HELL code by Bonvini for the small-x
*     resummation.
*
************************************************************************
      subroutine initHELL(la,pt,asmc,asmb,asmt)
*
      implicit none
*
      include "../commons/consts.h"
**
*     Input Variables
*
      integer la,pt
      double precision asmc,asmb,asmt
**
*     Internal Variables
*
      double precision alphasmc,alphasmb,alphasmt
*
      alphasmc = 4d0 * pi * asmc
      alphasmb = 4d0 * pi * asmb
      alphasmt = 4d0 * pi * asmt
*
      call HELLLogOrder(la)
      call HELLOrder(pt)
      call HELL(alphasmc,alphasmb,alphasmt)
*
      return
      end
