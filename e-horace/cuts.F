#ifndef REPORT
** generic cuts routine, to be modified and used by the user. Active only if "REPORT" is not defined.
      subroutine cuts(p3,p4,qph,icut)      
!! routine where experimental cuts are applied
* Inputs:
* ======
* p3(0:3) = (E,px,py,pz) = lepton momentum
* p4(0:3) = (E,px,py,pz) = neutrino momentum (lepton momentum for Z production)
* qph(40,0:3) = momenta of the emitted photons
*
* Outputs:
* =======
* icut = 0 ---> event accepted
* icut = 1 ---> event rejected
      include 'shared.inc'
      double precision inv_or_trans_mass_min
      common/fordebugging/idebugging
      common/yminimo/ymin ! for x1,x2pdf sampling
      common/xpdfminmax/xamin,xamax,xbmin,xbmax
      data xamin,xamax,xbmin,xbmax /1.d0,0.d0,1.d0,0.d0/
      dimension ptmp1(0:3),ptmp2(0:3),ptmp3(0:3)
      common/cutscommon/widercut,ifirst
      data ifirst,widercut /0,0.d0/

      common/select_partial_xs/ipartial1,ipartial2
      
      if (ifirst.eq.0) then
         
         pt1_min = 25d0  ! GeV
         pt2_min = 25d0  ! GeV
         eta_max = 2.5

!     invariant mass cut for g/Z production or transverse mass cut for W production
         inv_or_trans_mass_min = 15.d0 ! GeV

** do not change these lines [[
         ifirst    = 1
         ptmin     = pt1_min
         ptminmiss = pt2_min
         etamax    = eta_max
         tmmin     = inv_or_trans_mass_min
         ymin = max(((ptmin+ptminmiss)/2.d0/ebeam)**2,
     >               (tmmin/2.d0/ebeam)**2)
         ymin = max(ymin,xpdfmin*xpdfmin)
**                           ]]
      endif
      
      icut = 1 ! event rejected

** PT CUT *****************************************
      ptl    = sqrt(p3(1)**2+p3(2)**2)
      ptmiss = sqrt(p4(1)**2+p4(2)**2)
      if (ptl.lt.ptmin) then
         return
      endif
      if (ptmiss.lt.ptminmiss) then
         return
      endif
** RAPIDITY CUT ****************************************
      eta1 = abs(pseudorapidity(p3))
      eta2 = abs(pseudorapidity(p4))
      if (eta1.gt.etamax) then
         return
      endif
      if (boson.eq.'Z') then
         if (eta2.gt.etamax) return
      endif
*---- MT CUT for W
      if (tmmin.gt.0.d0.and.boson.eq.'W') then
         ptmp1(3) = 0.d0
         ptmp2(3) = 0.d0
         do k = 1,2
            ptmp1(k) = p3(k)
            ptmp2(k) = p4(k)
         enddo
         ptmp1(0) = sqrt(ptmp1(1)**2+ptmp1(2)**2)
         ptmp2(0) = sqrt(ptmp2(1)**2+ptmp2(2)**2)
         do k = 0,3
            ptmp3(k) = ptmp1(k) + ptmp2(k)
         enddo
         tm = sqrt(abs(dot(ptmp3,ptmp3)))
         if (tm.lt.tmmin) return
      endif
*---- INVARIANT MASS CUT for Z
      if (tmmin.gt.0.d0.and.boson.eq.'Z') then
         ptmp1 = p3 + p4
         if (sqrt(dot(ptmp1,ptmp1)).lt.tmmin) return
      endif

      icut = 0 ! event accepted
      return
      end
*********************************************************************************
#else
*********************************************************************************
** cut routine to reproduce report numbers, not to be modified by the user.
** Active only if "REPORT" is defined.
      subroutine cuts(p3,p4,qph,icut)      
      include 'shared.inc'
      double precision inv_or_trans_mass_min
      common/fordebugging/idebugging
      common/yminimo/ymin ! for x1,x2pdf sampling
      common/xpdfminmax/xamin,xamax,xbmin,xbmax
      data xamin,xamax,xbmin,xbmax /1.d0,0.d0,1.d0,0.d0/
      dimension ptmp1(0:3),ptmp2(0:3),ptmp3(0:3)
      common/cutscommon/widercut,ifirst
      data ifirst,widercut /0,0.d0/

      common/select_partial_xs/ipartial1,ipartial2
      
      if (ifirst.eq.0) then

         pt1_min = 25d0  ! GeV
         pt2_min = 25d0  ! GeV
         eta_max = 2.5

         inv_or_trans_mass_min = 0.d0 ! GeV

         if (boson.eq.'Z') then
            inv_or_trans_mass_min = 50.d0
         else  ! if W
            if (ischemealpha.eq.1) inv_or_trans_mass_min = 40.d0
         endif


** do not change these lines [[
         ifirst    = 1
         ptmin     = pt1_min
         ptminmiss = pt2_min
         etamax    = eta_max
         tmmin     = inv_or_trans_mass_min
         ymin = max(((ptmin+ptminmiss)/2.d0/ebeam)**2,
     >               (tmmin/2.d0/ebeam)**2)
         ymin = max(ymin,xpdfmin*xpdfmin)
**                           ]]
      endif
      
      icut = 1 ! event rejected

******   HARD WIRED CUT, for the report
      if (boson.eq.'W') then
         if (dot(p3+p4,p3+p4).lt.1.d0) return
      endif

******************************************
** PT CUT ********************************
      ptl    = sqrt(p3(1)**2+p3(2)**2)
      ptmiss = sqrt(p4(1)**2+p4(2)**2)
      if (ptl.lt.ptmin) then
         return
      endif
      if (ptmiss.lt.ptminmiss) then
         return
      endif
** RAPIDITY CUT **************************
      eta1 = abs(pseudorapidity(p3))
      eta2 = abs(pseudorapidity(p4))
      if (eta1.gt.etamax) then
         return
      endif
      if (boson.eq.'Z') then
         if (eta2.gt.etamax) return
      endif
*---- MT CUT for W
      if (tmmin.gt.0.d0.and.boson.eq.'W') then
         ptmp1(3) = 0.d0
         ptmp2(3) = 0.d0
         do k = 1,2
            ptmp1(k) = p3(k)
            ptmp2(k) = p4(k)
         enddo
         ptmp1(0) = sqrt(ptmp1(1)**2+ptmp1(2)**2)
         ptmp2(0) = sqrt(ptmp2(1)**2+ptmp2(2)**2)
         do k = 0,3
            ptmp3(k) = ptmp1(k) + ptmp2(k)
         enddo
         tm = sqrt(abs(dot(ptmp3,ptmp3)))
         if (tm.lt.tmmin) return
      endif
*---- INVARIANT MASS CUT for Z
      if (tmmin.gt.0.d0.and.boson.eq.'Z') then
         ptmp1 = p3 + p4
         if (sqrt(dot(ptmp1,ptmp1)).lt.tmmin) return
      endif

      icut = 0 ! event accepted
      return
      end
****************************************************
#endif

      function get_deltaR(p1,p2)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3)

! from ALPGEN
      pt1 = sqrt(p1(1)**2+p1(2)**2)
      pt2 = sqrt(p2(1)**2+p2(2)**2)
      eta1 = -log(tan(0.5d0*atan2(pt1,p1(3))))
      eta2 = -log(tan(0.5d0*atan2(pt2,p2(3))))
      dphi = (p1(1)*p2(1)+p1(2)*p2(2))/pt1/pt2
      if (abs(dphi).gt.1.d0) then
         dphi = dphi/abs(dphi)
      endif
      dphi = acos(dphi)

      get_deltaR = sqrt((eta1-eta2)**2+dphi**2)
      
      return
      end
**********

*********************
      subroutine recombine_e_phot(p3,p4,qph)
      include 'shared.inc'
      character*6 ord
      common/qedORDER/ord
      common/n_phot/nphot

      if (lepton.gt.1.or.ord.eq.'born'.or.nphot.eq.0) return

      do k = 1,nphot
         drp3 = get_deltaR(p3,qph(k,:))
         drp4 = get_deltaR(p4,qph(k,:))         
         icl = 3
         if (drp4.lt.drp3.and.boson.eq.'Z') icl = 4
         if (drp3.lt.0.1d0.and.icl.eq.3) then
            p3 = p3 + qph(k,:)
            qph(k,:) = 0.d0
         endif
         if (drp4.lt.0.1d0.and.icl.eq.4) then
            p4 = p4 + qph(k,:)
            qph(k,:) = 0.d0
         endif
      enddo
      
      return
      end
