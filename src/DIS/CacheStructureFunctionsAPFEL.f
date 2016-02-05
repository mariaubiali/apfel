************************************************************************
*
*     CacheStructureFunctionsAPFEL.f
*
*     This subroutine chaces the evolved PDFs over a grid in x and Q2
*     in such a way that the DGLAP evolution is done offline only once
*     and PDFs at generic values of Q are obtained by interpolation.
*
************************************************************************
      subroutine CacheStructureFunctionsAPFEL(Q0)
*
      implicit none
*
      include "../commons/Welcome.h"
      include "../commons/grid.h"
      include "../commons/gridQ.h"
      include "../commons/Evs.h"
      include "../commons/PDFEvolution.h"
      include "../commons/scales.h"
      include "../commons/MaxFlavourPDFs.h"
      include "../commons/MaxFlavourAlpha.h"
      include "../commons/pdfset.h"
      include "../commons/StructureFunctions.h"
      include "../commons/StructureFunctionsxQ.h"
      include "../commons/m2th.h"
      include "../commons/Nf_FF.h"
      include "../commons/TMC.h"
**
*     Input Variables
*
      double precision Q0
**
*     Internal Variables
*
      integer nfmax
      integer isg
      integer iq2,iq2c
      integer alpha
      integer iproc
      integer icomp
      double precision lnQmin,lnQmax
      double precision eps
      double precision t1,t2
      character*100 pdfsetbkp
      parameter(eps=1d-14)
*
      call cpu_time(t1)
*
*     Define maximun number of flavours
*
      nfmax = max(nfMaxPDFs,nfMaxAlpha)
*
*     Compute Q2 grid.
*     Use a distribution of the Q2 nodes uniform in ln(ln(Q2/Lam2))
*     and that has nodes on the heavy quark thresholds.
*
      lnQmin = dlog( Q2min / Lam2 )
      lnQmax = dlog( Q2max / Lam2 )
*
*     Initialize number of points per subgrid
*
      do isg=3,7
         nQ(isg) = 0
      enddo
*
      if(Evs.eq."VF")then
         if(Q2min.ge.m2th(6))then
            nfin = 6
         elseif(Q2min.ge.m2th(5))then
            nfin = 5
         elseif(Q2min.ge.m2th(4))then
            nfin = 4
         else
            nfin = 3
         endif
         if(nfin.gt.nfmax) nfin = nfmax
*
         if(Q2max.gt.m2th(6))then
            nffi = 6
         elseif(Q2max.gt.m2th(5))then
            nffi = 5
         elseif(Q2max.gt.m2th(4))then
            nffi = 4
         else
            nffi = 3
         endif
         if(nffi.gt.nfmax) nffi = nfmax
*
         isg = nfin
         do iq2=0,nQ2g
            Q2g(iq2) = Lam2 * dexp( lnQmin
     1                 * dexp( dble( iq2 ) / dble( nQ2g )
     2                 * dlog( lnQmax / lnQmin ) ) )
            if(Q2g(iq2).lt.m2th(isg+1)+eps)then
               nQ(isg) = nQ(isg) + 1
            else
               isg = isg + 1
               if(isg.gt.nfmax) isg = nfmax
               nQ(isg) = nQ(isg) + 1
            endif
         enddo
*
*     Make sure that all subgrids have at least two points.
*
         do isg=nfin,nffi-1
            if(nQ(isg).lt.2)then
               nQ(isg) = nQ(isg) + 1
               nQ(isg+1) = nQ(isg+1) - 1
            endif
         enddo
*
*     Redefine the grid with the subgrids
*
         lnQmin = dlog( Q2min / Lam2 )
         if(nfin.eq.nffi)then
            lnQmax = dlog( Q2max / Lam2 )
         else
            lnQmax = dlog( m2th(nfin+1) / Lam2 )
         endif
*
         iq2c = -1
         do isg=nfin,nffi
            do iq2=1,nQ(isg)
               iq2c = iq2c + 1
               Q2g(iq2c) = Lam2 * dexp( lnQmin
     1              * dexp( dble( iq2 - 1 )
     2              / dble( nQ(isg) - 1 )
     3              * dlog( lnQmax / lnQmin ) ) )
            enddo
            lnQmin = dlog( ( 1d0 + eps ) * Q2g(iq2c) / Lam2 )
            if(isg.eq.nffi-1)then
               lnQmax = dlog( Q2max / Lam2 )
            else
               lnQmax = dlog( m2th(isg+2) / Lam2 )
            endif
            Q2g(iq2c) = ( 1d0 - eps ) * Q2g(iq2c)
         enddo
*
         if(iq2c.ne.nQ2g)then
            write(6,*) "In CacheStructureFunctionsAPFEL.f:"
            write(6,*) "Mismatch in the Number of Q2 nodes"
            write(6,*) "- Expected = ",nQ2g+1
            write(6,*) "- Found = ",iq2c+1
            call exit(-10)
         endif
      else
         nfin = Nf_FF
         nffi = Nf_FF
*
         nQ(Nf_FF) = nQ2g
*
         do iq2=0,nQ2g
            Q2g(iq2) = Lam2 * dexp( lnQmin
     1                 * dexp( dble( iq2 ) / dble( nQ2g )
     2                 * dlog( lnQmax / lnQmin ) ) )
         enddo
      endif
*     Backup PDF name
      pdfsetbkp = pdfset
*
*     Evolve PDFs on the Q2-grid
*
      Q2g(-1) = Q0 * Q0
      do iq2=0,nQ2g
         do iproc=0,2
            if(iproc.eq.0) call SetProcessDIS("EM")
            if(iproc.eq.1) call SetProcessDIS("NC")
            if(iproc.eq.2) call SetProcessDIS("CC")
            if(Q0.lt.0d0.or.iproc.gt.0)then
               call ComputeStructureFunctionsAPFEL(dsqrt(Q2g(iq2)),
     1                                             dsqrt(Q2g(iq2)))
            else
               if(PDFevol(1:9).eq."truncated")then
                  call ComputeStructureFunctionsAPFEL(Q0,
     1                                                dsqrt(Q2g(iq2)))
               else
                  call ComputeStructureFunctionsAPFEL(dsqrt(Q2g(iq2-1)),
     1                                                dsqrt(Q2g(iq2)))
                  call SetPDFSet("apfel")
               endif
            endif
            if(TMC)then
               
            else
               do alpha=0,nin(0)
                  do icomp=3,7
                     SFxQ(iproc,1,icomp,alpha,iq2) = F2(icomp,0,alpha)
                     SFxQ(iproc,2,icomp,alpha,iq2) = FL(icomp,0,alpha)
                     SFxQ(iproc,3,icomp,alpha,iq2) = F3(icomp,0,alpha)
                  enddo
c                  write(6,*) "... ",iproc,dsqrt(Q2g(iq2),
c     1                 SFxQ(iproc,1,6,alpha,iq2)
               enddo
            endif
         enddo
      enddo
*     Restore PDF name and reset PDFs at the initial scale
      call SetPDFSet(pdfsetbkp)
      call ComputeStructureFunctionsAPFEL(Q0,Q0)
*
*     Caching complete
*
      InCacheSFs = "done"
*
      call cpu_time(t2)
*
      if(Welcome)then
         write(6,*) achar(27)//"[34m"//"Structure functions",
     1                         " have been cached"
         write(6,*) achar(27)//"[0m"
         write(6,"(a,f7.3,a)") " Caching completed in",t2-t1," s"
         write(6,*) " "
      endif
*
      return
      end