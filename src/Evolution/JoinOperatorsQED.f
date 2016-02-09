************************************************************************
*
*     JoinOperatorsQED.f:
*
*     This routine joins the QEDxQCD evolution operators computed with 
*     different numbers of active flavours.
*
*     QEDxQCD evolution basis:
*     0   1   2   3   4   5   6   7   8   9  10  11  12  13
*     gm  Sg   g   V  V3  V8 V15 V24 V35  T3  T8 T15 T24 T35
*     
************************************************************************
      subroutine JoinOperatorsQED(jgrid)
*
      implicit none
*
      include "../commons/grid.h"
      include "../commons/ThresholdAlphaQCD.h"
      include "../commons/EvolutionMatrices.h"
      include "../commons/transQCD.h"
      include "../commons/EvolutionOperator.h"
**
*     Input Variables
*
      integer jgrid
**
*     Internal Variables
*
      integer i,j,k,l
      integer nf,nfm
      integer alpha,beta,gamma,delta
      double precision coup,integralsMatching
      double precision MatQCDns(0:nint_max,0:nint_max),Match(2)
      double precision MatQCDsg(2,2,0:nint_max,0:nint_max)
      double precision EvQCDb(0:13,0:13,0:nint_max,0:nint_max)
      double precision EvQCD(0:13,0:13,0:nint_max,0:nint_max)
*
      if(sgn.ne.1)then
         write(6,*) "In JoinOperatorsQED.f:"
         write(6,*) "Backward evolution not allowed yet."
         call exit(-10)
      endif
*
      do alpha=0,nin(igrid)
         do beta=alpha,nin(igrid)
*     Set evolution operators to zero
            do i=0,13
               do j=0,13
                  EvQCD(i,j,alpha,beta) = 0d0
               enddo
            enddo
*     Singlet 1 (excluding leptons)
            do i=1,4
               do j=1,4
                  EvQCD(i-1,j-1,alpha,beta) =
     1                 MUnisg1(nfi,nli,i,j,alpha,beta)
               enddo
            enddo
*     Singlet 2
            do i=1,2
               do j=1,2
                  EvQCD(i+3,j+3,alpha,beta) =
     1                 MUnisg2(nfi,nli,i,j,alpha,beta)
               enddo
            enddo
*     Tu1
            if(nfi.ge.4)then
               EvQCD(6,6,alpha,beta) =
     1              MUninspu(nfi,nli,alpha,beta)
            else
               do i=1,5
                  EvQCD(6,6,alpha,beta) =
     1                 ( MUnisg1(nfi,nli,3,i,alpha,beta)
     2                 + MUnisg1(nfi,nli,4,i,alpha,beta) ) / 2d0
               enddo
            endif
*     Tu2
            if (nfi.ge.6)then
               EvQCD(7,7,alpha,beta) =
     1              MUninspd(nfi,nli,alpha,beta)
            else
               do i=1,5
                  EvQCD(7,7,alpha,beta) =
     1                 ( MUnisg1(nfi,nli,3,i,alpha,beta)
     2                 + MUnisg1(nfi,nli,4,i,alpha,beta) ) / 2d0
               enddo
            endif
*     Td1
            EvQCD(8,8,alpha,beta) =
     1           MUninspd(nfi,nli,alpha,beta)
*     Td2
            if(nfi.ge.5)then
               EvQCD(9,9,alpha,beta) =
     1              MUninspd(nfi,nli,alpha,beta)
            else
               do i=1,5
                  EvQCD(9,9,alpha,beta) =
     1                 ( MUnisg1(nfi,nli,3,i,alpha,beta)
     2                 - MUnisg1(nfi,nli,4,i,alpha,beta) ) / 2d0
               enddo
            endif
*     Vu1
            if(nfi.ge.4)then
               EvQCD(10,10,alpha,beta) =
     1              MUninsmu(nfi,nli,alpha,beta)
            else
               do i=1,2
                  EvQCD(10,10,alpha,beta) =
     1                 ( MUnisg2(nfi,nli,1,i,alpha,beta)
     2                 + MUnisg2(nfi,nli,2,i,alpha,beta) ) / 2d0
               enddo               
            endif
*     Vu2
            if(nfi.ge.6)then
               EvQCD(11,11,alpha,beta) =
     1              MUninsmu(nfi,nli,alpha,beta)
            else
               do i=1,2
                  EvQCD(11,11,alpha,beta) =
     1                 ( MUnisg2(nfi,nli,1,i,alpha,beta)
     2                 + MUnisg2(nfi,nli,2,i,alpha,beta) ) / 2d0
               enddo
            endif
*     Vd1
            EvQCD(12,12,alpha,beta) =
     1           MUninsmd(nfi,nli,alpha,beta)
*     Vd2
            if(nfi.ge.5)then
               EvQCD(13,13,alpha,beta) =
     1              MUninsmd(nfi,nli,alpha,beta)
            else
               do i=1,2
                  EvQCD(13,13,alpha,beta) =
     1                 ( MUnisg2(nfi,nli,1,i,alpha,beta)
     2                 - MUnisg2(nfi,nli,2,i,alpha,beta) ) / 2d0
               enddo
            endif            
         enddo
      enddo
c$$$*
c$$$*     If the initial and the final numbers of flavours are different ...
c$$$*
c$$$      if(nfi.ne.nff)then
c$$$         do nf=nfi,nff-1
c$$$*
c$$$*     Set temporary evolution operators to zero
c$$$*
c$$$            do alpha=0,nin(igrid)
c$$$               do beta=alpha,nin(igrid)
c$$$                  do i=0,13
c$$$                     do j=0,13
c$$$                        EvQCDb(i,j,alpha,beta) = 0d0
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$*     
c$$$            nfm = nf + 1
c$$$*     Get alphas value at the heavy quark threshold (with nfm active flavours)
c$$$c            coup = asthUp(nfm)
c$$$            coup = asthDown(nfm)
c$$$*     Contruct matching conditions at this threshod
c$$$            if(IsExt(igrid))then
c$$$               do alpha=0,nin(igrid)
c$$$                  do beta=alpha,nin(igrid)
c$$$                     MatQCDns(alpha,beta)     =
c$$$     1                    integralsMatching(nf+1,alpha,beta,coup,1,sgn)
c$$$                     MatQCDsg(1,1,alpha,beta) =
c$$$     1                    integralsMatching(nf+1,alpha,beta,coup,2,sgn)
c$$$                     MatQCDsg(1,2,alpha,beta) =
c$$$     1                    integralsMatching(nf+1,alpha,beta,coup,3,sgn)
c$$$                     MatQCDsg(2,1,alpha,beta) =
c$$$     1                    integralsMatching(nf+1,alpha,beta,coup,4,sgn)
c$$$                     MatQCDsg(2,2,alpha,beta) =
c$$$     1                    integralsMatching(nf+1,alpha,beta,coup,5,sgn)
c$$$                  enddo
c$$$               enddo
c$$$            else
c$$$               do alpha=0,nin(igrid)
c$$$                  do beta=alpha,nin(igrid)
c$$$                     MatQCDns(alpha,beta)     =
c$$$     1                   integralsMatching(nf+1,0,beta-alpha,coup,1,sgn)
c$$$                     MatQCDsg(1,1,alpha,beta) =
c$$$     1                   integralsMatching(nf+1,0,beta-alpha,coup,2,sgn)
c$$$                     MatQCDsg(1,2,alpha,beta) =
c$$$     1                   integralsMatching(nf+1,0,beta-alpha,coup,3,sgn)
c$$$                     MatQCDsg(2,1,alpha,beta) =
c$$$     1                   integralsMatching(nf+1,0,beta-alpha,coup,4,sgn)
c$$$                     MatQCDsg(2,2,alpha,beta) =
c$$$     1                   integralsMatching(nf+1,0,beta-alpha,coup,5,sgn)
c$$$                  enddo
c$$$               enddo
c$$$            endif
c$$$*     
c$$$*     Now combine evolution tables for different numbers of active
c$$$*     flavours.
c$$$*     
c$$$            do alpha=0,nin(igrid)
c$$$               do beta=alpha,nin(igrid)
c$$$*     Singlet and Gluon
c$$$                  do i=1,2
c$$$                     do j=1,2
c$$$                        do gamma=0,nin(igrid)
c$$$                           do delta=0,nin(igrid)
c$$$                              do k=1,2
c$$$                                 do l=1,2
c$$$                                    EvQCDb(i,j,alpha,beta) = 
c$$$     1                                   EvQCDb(i,j,alpha,beta)
c$$$     2                                   + MQCDsg(nf+1,i,k,alpha,gamma)
c$$$     3                                   * MatQCDsg(k,l,gamma,delta)
c$$$     4                                   * EvQCD(l,j,delta,beta)
c$$$                                 enddo
c$$$                              enddo
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$                  enddo
c$$$*     
c$$$                  do gamma=0,nin(igrid)
c$$$                     do delta=0,nin(igrid)
c$$$*     Total Valence
c$$$                        EvQCDb(3,3,alpha,beta) = 
c$$$     1                       EvQCDb(3,3,alpha,beta)
c$$$     2                       + MQCDnsv(nf+1,alpha,gamma)
c$$$     3                       * MatQCDns(gamma,delta)
c$$$     4                       * EvQCD(3,3,delta,beta)
c$$$*     V3
c$$$                        EvQCDb(4,4,alpha,beta) = 
c$$$     1                       EvQCDb(4,4,alpha,beta)
c$$$     2                       + MQCDnsm(nf+1,alpha,gamma)
c$$$     3                       * MatQCDns(gamma,delta)
c$$$     4                       * EvQCD(4,4,delta,beta)
c$$$*     V8
c$$$                        EvQCDb(5,5,alpha,beta) = 
c$$$     1                       EvQCDb(5,5,alpha,beta)
c$$$     2                       + MQCDnsm(nf+1,alpha,gamma)
c$$$     3                       * MatQCDns(gamma,delta)
c$$$     4                       * EvQCD(5,5,delta,beta)
c$$$*     T3
c$$$                        EvQCDb(9,9,alpha,beta) = 
c$$$     1                       EvQCDb(9,9,alpha,beta)
c$$$     2                       + MQCDnsp(nf+1,alpha,gamma)
c$$$     3                       * MatQCDns(gamma,delta)
c$$$     4                       * EvQCD(9,9,delta,beta)
c$$$*     T8
c$$$                        EvQCDb(10,10,alpha,beta) = 
c$$$     1                       EvQCDb(10,10,alpha,beta)
c$$$     2                       + MQCDnsp(nf+1,alpha,gamma)
c$$$     3                       * MatQCDns(gamma,delta)
c$$$     4                       * EvQCD(10,10,delta,beta)
c$$$                     enddo
c$$$                  enddo
c$$$*     Charm threshold
c$$$                  if(nfm.eq.4)then
c$$$                     do gamma=0,nin(igrid)
c$$$                        do delta=0,nin(igrid)
c$$$*     V15
c$$$                           EvQCDb(6,3,alpha,beta) = 
c$$$     1                          EvQCDb(6,3,alpha,beta)
c$$$     2                          + MQCDnsm(nf+1,alpha,gamma)
c$$$     3                          * MatQCDns(gamma,delta)
c$$$     4                          * EvQCD(6,3,delta,beta)
c$$$*     V24
c$$$                           EvQCDb(7,3,alpha,beta) = 
c$$$     1                          EvQCDb(7,3,alpha,beta)
c$$$     2                          + MQCDnsv(nf+1,alpha,gamma)
c$$$     3                          * MatQCDns(gamma,delta)
c$$$     4                          * EvQCD(7,3,delta,beta)
c$$$*     V35
c$$$                           EvQCDb(8,3,alpha,beta) = 
c$$$     1                          EvQCDb(8,3,alpha,beta)
c$$$     2                          + MQCDnsv(nf+1,alpha,gamma)
c$$$     3                          * MatQCDns(gamma,delta)
c$$$     4                          * EvQCD(8,3,delta,beta)
c$$$*     T15
c$$$                           Match(1) = MatQCDns(gamma,delta) 
c$$$     1                          - 3d0 * ( MatQCDsg(1,1,gamma,delta) 
c$$$     2                          - MatQCDns(gamma,delta) )
c$$$                           Match(2) = - 3d0 * MatQCDsg(1,2,gamma,delta)
c$$$                           do j=1,2
c$$$                              do k=1,2
c$$$                                 EvQCDb(11,j,alpha,beta) = 
c$$$     1                                EvQCDb(11,j,alpha,beta)
c$$$     2                                + MQCDnsp(nf+sgn,alpha,gamma)
c$$$     3                                * Match(k)
c$$$     4                                * EvQCD(k,j,delta,beta)
c$$$                              enddo
c$$$                           enddo
c$$$*     
c$$$                           do j=1,2
c$$$*     T24
c$$$                              do k=1,2
c$$$                                 do l=1,2
c$$$                                    EvQCDb(12,j,alpha,beta) = 
c$$$     1                                   EvQCDb(12,j,alpha,beta)
c$$$     2                                   + MQCDsg(nf+1,1,k,alpha,gamma)
c$$$     3                                   * MatQCDsg(k,l,gamma,delta)
c$$$     4                                   * EvQCD(l,j,delta,beta)
c$$$                                 enddo
c$$$                              enddo
c$$$*     T35
c$$$                              EvQCDb(13,j,alpha,beta) = 
c$$$     1                             EvQCDb(12,j,alpha,beta)
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$*     Bottom threshold
c$$$                  elseif(nfm.eq.5)then
c$$$                     do gamma=0,nin(igrid)
c$$$                        do delta=0,nin(igrid)
c$$$*     V15
c$$$                           if(nfi.ge.4)then
c$$$                              EvQCDb(6,6,alpha,beta) = 
c$$$     1                             EvQCDb(6,6,alpha,beta)
c$$$     2                             + MQCDnsm(nf+1,alpha,gamma)
c$$$     3                             * MatQCDns(gamma,delta)
c$$$     4                             * EvQCD(6,6,delta,beta)
c$$$                           else
c$$$                              EvQCDb(6,3,alpha,beta) = 
c$$$     1                             EvQCDb(6,3,alpha,beta)
c$$$     2                             + MQCDnsm(nf+1,alpha,gamma)
c$$$     3                             * MatQCDns(gamma,delta)
c$$$     4                             * EvQCD(6,3,delta,beta)
c$$$                           endif
c$$$*     V24
c$$$                           EvQCDb(7,3,alpha,beta) = 
c$$$     1                          EvQCDb(7,3,alpha,beta)
c$$$     2                          + MQCDnsm(nf+1,alpha,gamma)
c$$$     3                          * MatQCDns(gamma,delta)
c$$$     4                          * EvQCD(7,3,delta,beta)
c$$$*     V35
c$$$                           EvQCDb(8,3,alpha,beta) = 
c$$$     1                          EvQCDb(8,3,alpha,beta)
c$$$     2                          + MQCDnsv(nf+1,alpha,gamma)
c$$$     3                          * MatQCDns(gamma,delta)
c$$$     4                          * EvQCD(8,3,delta,beta)
c$$$*     T15
c$$$                           if(nfi.ge.4)then
c$$$                              EvQCDb(11,11,alpha,beta) = 
c$$$     1                             EvQCDb(11,11,alpha,beta)
c$$$     2                             + MQCDnsp(nf+1,alpha,gamma)
c$$$     3                             * MatQCDns(gamma,delta)
c$$$     4                             * EvQCD(11,11,delta,beta)
c$$$                           else
c$$$                              do j=1,2
c$$$                                 EvQCDb(11,j,alpha,beta) = 
c$$$     1                                EvQCDb(11,j,alpha,beta)
c$$$     2                                + MQCDnsp(nf+1,alpha,gamma)
c$$$     3                                * MatQCDns(gamma,delta)
c$$$     4                                * EvQCD(11,j,delta,beta)
c$$$                              enddo
c$$$                           endif
c$$$*     T24
c$$$                           Match(1) = MatQCDns(gamma,delta) 
c$$$     1                          - 4d0 * ( MatQCDsg(1,1,gamma,delta) 
c$$$     2                          - MatQCDns(gamma,delta) )
c$$$                           Match(2) = - 4d0 * MatQCDsg(1,2,gamma,delta)
c$$$                           do j=1,2
c$$$                              do k=1,2
c$$$                                 EvQCDb(12,j,alpha,beta) = 
c$$$     1                                EvQCDb(12,j,alpha,beta)
c$$$     2                                + MQCDnsp(nf+1,alpha,gamma)
c$$$     3                                * Match(k)
c$$$     4                                * EvQCD(k,j,delta,beta)
c$$$                              enddo
c$$$                           enddo
c$$$*     T35
c$$$                           do j=1,2
c$$$                              do k=1,2
c$$$                                 do l=1,2
c$$$                                    EvQCDb(13,j,alpha,beta) = 
c$$$     1                                   EvQCDb(13,j,alpha,beta)
c$$$     2                                   + MQCDsg(nf+1,1,k,alpha,gamma)
c$$$     3                                   * MatQCDsg(k,l,gamma,delta)
c$$$     4                                   * EvQCD(l,j,delta,beta)
c$$$                                 enddo
c$$$                              enddo
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$*     Top threshold
c$$$                  elseif(nfm.eq.6)then
c$$$                     do gamma=0,nin(igrid)
c$$$                        do delta=0,nin(igrid)
c$$$*     V15
c$$$                           if(nfi.ge.4)then
c$$$                              EvQCDb(6,6,alpha,beta) = 
c$$$     1                             EvQCDb(6,6,alpha,beta)
c$$$     2                             + MQCDnsm(nf+1,alpha,gamma)
c$$$     3                             * MatQCDns(gamma,delta)
c$$$     4                             * EvQCD(6,6,delta,beta)
c$$$                           else
c$$$                              EvQCDb(6,3,alpha,beta) = 
c$$$     1                             EvQCDb(6,3,alpha,beta)
c$$$     2                             + MQCDnsm(nf+1,alpha,gamma)
c$$$     3                             * MatQCDns(gamma,delta)
c$$$     4                             * EvQCD(6,3,delta,beta)
c$$$                           endif
c$$$*     V24
c$$$                           if(nfi.ge.5)then
c$$$                              EvQCDb(7,7,alpha,beta) = 
c$$$     1                             EvQCDb(7,7,alpha,beta)
c$$$     2                             + MQCDnsm(nf+1,alpha,gamma)
c$$$     3                             * MatQCDns(gamma,delta)
c$$$     4                             * EvQCD(7,7,delta,beta)
c$$$                           else
c$$$                              EvQCDb(7,3,alpha,beta) = 
c$$$     1                             EvQCDb(7,3,alpha,beta)
c$$$     2                             + MQCDnsm(nf+1,alpha,gamma)
c$$$     3                             * MatQCDns(gamma,delta)
c$$$     4                             * EvQCD(7,3,delta,beta)
c$$$                           endif
c$$$*     V35
c$$$                           EvQCDb(8,3,alpha,beta) = 
c$$$     1                          EvQCDb(8,3,alpha,beta)
c$$$     2                          + MQCDnsm(nf+1,alpha,gamma)
c$$$     3                          * MatQCDns(gamma,delta)
c$$$     4                          * EvQCD(8,3,delta,beta)
c$$$*     T15
c$$$                           if(nfi.ge.4)then
c$$$                              EvQCDb(11,11,alpha,beta) = 
c$$$     1                             EvQCDb(11,11,alpha,beta)
c$$$     2                             + MQCDnsp(nf+1,alpha,gamma)
c$$$     3                             * MatQCDns(gamma,delta)
c$$$     4                             * EvQCD(11,11,delta,beta)
c$$$                           else
c$$$                              do j=1,2
c$$$                                 EvQCDb(11,j,alpha,beta) = 
c$$$     1                                EvQCDb(11,j,alpha,beta)
c$$$     2                                + MQCDnsp(nf+1,alpha,gamma)
c$$$     3                                * MatQCDns(gamma,delta)
c$$$     4                                * EvQCD(11,j,delta,beta)
c$$$                              enddo
c$$$                           endif
c$$$*     T24
c$$$                           if(nfi.ge.5)then
c$$$                              EvQCDb(12,12,alpha,beta) = 
c$$$     1                             EvQCDb(12,12,alpha,beta)
c$$$     2                             + MQCDnsp(nf+1,alpha,gamma)
c$$$     3                             * MatQCDns(gamma,delta)
c$$$     4                             * EvQCD(12,12,delta,beta)
c$$$                           else
c$$$                              do j=1,2
c$$$                                 EvQCDb(12,j,alpha,beta) = 
c$$$     1                                EvQCDb(12,j,alpha,beta)
c$$$     2                                + MQCDnsp(nf+1,alpha,gamma)
c$$$     3                                * MatQCDns(gamma,delta)
c$$$     4                                * EvQCD(12,j,delta,beta)
c$$$                              enddo
c$$$                           endif
c$$$*     T35
c$$$                           Match(1) = MatQCDns(gamma,delta) 
c$$$     1                          - 5d0 * ( MatQCDsg(1,1,gamma,delta) 
c$$$     2                          - MatQCDns(gamma,delta) )
c$$$                           Match(2) = - 5d0 * MatQCDsg(1,2,gamma,delta)
c$$$                           do j=1,2
c$$$                              do k=1,2
c$$$                                 EvQCDb(13,j,alpha,beta) = 
c$$$     1                                EvQCDb(13,j,alpha,beta)
c$$$     2                                + MQCDnsp(nf+1,alpha,gamma)
c$$$     3                                * Match(k)
c$$$     4                                * EvQCD(k,j,delta,beta)
c$$$                              enddo
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$                  endif
c$$$               enddo
c$$$            enddo
c$$$*     
c$$$*     Copy the backup evolution operators into the main ones
c$$$*     
c$$$            do alpha=0,nin(igrid)
c$$$               do beta=alpha,nin(igrid)
c$$$                  do i=0,13
c$$$                     do j=0,13
c$$$                        EvQCD(i,j,alpha,beta) = EvQCDb(i,j,alpha,beta)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      endif
c$$$*
c$$$*     Tranform the evolution operators from the evolution to the physical basis
c$$$*
c$$$      do alpha=0,nin(igrid)
c$$$         do beta=0,nin(igrid)
c$$$            do i=0,13
c$$$               do j=0,13
c$$$                  Ev2EvQCD(jgrid,i,j,alpha,beta) = EvQCD(i,j,alpha,beta)
c$$$                  Ev2PhQCD(jgrid,i-7,j,alpha,beta) = 0d0
c$$$                  Ph2PhQCD(jgrid,i-7,j-7,alpha,beta) = 0d0
c$$$                  do k=0,13
c$$$                     Ev2PhQCD(jgrid,i-7,j,alpha,beta) = 
c$$$     1                    Ev2PhQCD(jgrid,i-7,j,alpha,beta)
c$$$     2                       + Tev2phQCD(nff,i,k)
c$$$     3                       * EvQCD(k,j,alpha,beta)
c$$$                     do l=0,13
c$$$                        Ph2PhQCD(jgrid,i-7,j-7,alpha,beta) = 
c$$$     1                       Ph2PhQCD(jgrid,i-7,j-7,alpha,beta)
c$$$     2                       + Tev2phQCD(nff,i,k)
c$$$     3                       * EvQCD(k,l,alpha,beta)
c$$$     4                       * Tph2evQCD(nfi,l,j)
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      enddo

      return
      end
