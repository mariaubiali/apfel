************************************************************************
*
*     F2total.f:
*
*     This function returns the value of the inclusive structure function
*     F2.
*
************************************************************************
      function F2total(x)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/StructureFunctions.h"
**
*     Input Variables
*
      double precision x
**
*     Internal Variables
*
      integer n
      integer alpha
      double precision w_int
**
*     Output Variables
*
      double precision F2total
*
      if(x.lt.xmin(1).or.x.gt.xmax)then
         write(6,*) "In F2total.f:"
         write(6,*) "Invalid value of x =",x
         call exit(-10)
      endif
*
*     Select the grid
*
      do igrid=1,ngrid
         if(x.ge.xmin(igrid).and.x.lt.xmin(igrid+1))then
            goto 101
         endif
      enddo
*
*     Interpolation
*
 101  F2total = 0d0
      n = inter_degree(igrid)
      do alpha=0,nin(igrid)
         F2total = F2total + w_int(n,alpha,x) * F2(7,igrid,alpha)
      enddo
      if(dabs(F2total).le.1d-14) F2total = 0d0
*
      return
      end