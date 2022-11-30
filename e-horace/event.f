      subroutine generate_event_matched(sdif,i_ok,p3,p4,qph)
      include 'shared.inc'
      dimension pin1(0:3),pin2(0:3),qphcm(40,0:3)
      ! p1 --> u, p2 --> dbar, p3 --> neutrino, p4 --> positron
      double precision maaa,mbbb
      complex*16 mwcomplex
      complex*16 mzcomplex
      dimension vboost(0:3),q(0:3)
      dimension ama(4)
      dimension p3cm(0:3),p4cm(0:3),qcm(0:3),p1cm(0:3),p2cm(0:3),
     .     p5(0:3),p5cm(0:3)
      real*4 csi(1)
      character*10 model
      character*6 ord
      common/qedORDER/ord
      common/nphot_mode/nphotmode
      common/input_matching/model
      common/momentainitial/p1cm,p2cm
      common/masseinitial/maaa,mbbb
      common/eventinCOM/p3cm,p4cm,qphcm
      common/resonant/mwcomplex
      common /regulators/ESOFTMAX
      common/partialweights/wlepang,wphoten,wphotang,reg
      common/wparameters/amw,agw
      common/zparameters/amz,agz
      common/ifirstime_matched/ifirst_exewbck,ifirst_matched,ifirst_in,
     .     ifirst2
      data ifirst_matched /0/
      data ifirst_in,ifirst2 /0,0/
      common/chargesmasses_cm/charges(4),amasses(4)
      common/n_phot/nphot
      common/nrandom/ir
      common/sampling_weight/wsampling,spaziofasi,flux
      common/fordebugging/idebugging
      common/svfordebug/svdebug,emtx,emtxsub,wnphot
      common/cosinesforsubtraction/cosu,cosd
      common/scalapersottrazione/sqrtscalesub
      common/deltaerre/deltar
      common/rescaletogmu/dpu
      common/cosinesforcuts/cosenomax,cosenomin
      common/counter_event/icountevent
      data icountevent /0/
      common/subtractedsv/subtractedoalsv
      character*5 zIS
      common/selectvirtualz/zIS
      common/orderalphaphind/ioalphind
      common/ipartsmirrrored/ipart1mirr,ipart2mirr
      common/mcsamplingweight/w
      common/feyhorflags/dalpha,dalphaan,drho,drhoir1,drhoir12,rhofiu,
     .     rhofid,dkl,dku,dkd,ischemepippo
      common/forborncrosssection/phsp2b,flux2b,bornme,bornmeclean
      common/debugging/xxxx1,xxxx2

      common/select_partial_xs/ipartial1,ipartial2

      logical isnan
      
      icountevent = icountevent + 1
      sdif = 0.d0
      w         = 1.d0
      wsampling = 1.d0
      in_conf_spin = 4

ccc
      isubtract = iplsub
      if (ifirst_in.eq.0) then
         dpu   = alpha_gf*alpha_gf/alpha/alpha
         if (boson.eq.'Z'.and.ischemealpha.ne.1) dpu = 1.d0

         nphot = 1
c         call init_alpha
         if (isubtract.eq.1.and.i_pdf.eq.1.and.
     .        (ord.eq.'alpha'.or.ord.eq.'exp')) then
            if (boson.eq.'W') call init_deltaq(mw*rescalepdfscale)
            if (boson.eq.'Z') call init_deltaq(mz*rescalepdfscale)
         endif

         dalpha   = 0.d0
         dalphaan = 0.d0
         drhoir1  = 0.d0
         drhoir12 = 0.d0
         rhofiu   = 1.d0
         rhofid   = 1.d0
         dkl      = 0.d0
         dku      = 0.d0
         dkd      = 0.d0
         
         nphotmax = 0
         do i = 1,20
            sumsez(i)   = 0.d0
            sumsez2(i)  = 0.d0
            sezioni(i)  = 0.d0
            sezioni2(i) = 0.d0
            icounterpartial(i) = 0
         enddo
         mwcomplex = dcmplx(mw,-gw/2.d0)
         if (boson.eq.'Z') mwcomplex = dcmplx(mz,-gz/2.d0)
         amw = mw
         agw = gw
         amz = mz
         agz = gz
         subtractionpl = 0.d0
         ifirst_in = 1
      endif

! nullifying qph(1...nphot) as many as there were in the previous event...
      do j=0,3
         do k = 1,max(nphot,1)  ! nphot initialized to 1...            
            qph(k,j) = 0.d0
         enddo
         vboost(j) = 0.d0
         p5(j) = 0.d0
      enddo
      i_ok = 1

************* W+ W- production at the same time!!
      wwpwm = 1.d0
      if (abs(chfs).gt.0.d0.and.ichfs.eq.2) then
         chfs = 1.d0
         call wraprng(csi,1)
         if (csi(1).le.0.5) chfs = -1.d0
         wwpwm = 2.d0
      endif
**************************************
*********************************
      s = 4.d0*ebeam1*ebeam1
      call pdfsampling(s,spdf,wpdfpdfpdfjac)
      if (ipart1.eq.(ifl+1)) then
         i_ok = 0
         return
      endif

c      ipartial1 = 0
c      ipartial2 = 0
c      if (abs(ipartial1).lt.100.and.abs(ipartial2).lt.100) then
c         if (((ipart1.eq.ipartial1).and.(ipart2.eq.ipartial2)).or.
c     >       ((ipart1.eq.ipartial2).and.(ipart2.eq.ipartial1))) then
c            continue
c         else
c            i_ok = 0
c            return
c         endif
c      endif

      shat = s*x1pdf*x2pdf
      spdfcomm = spdf
*********************************

      colfac = 1.d0/3.d0
      if (ipart1.eq.0.and.ipart2.eq.0) colfac = 1.d0
** Notice: I associate Vud2 to pdf weight to enhance efficiency
      Vud2 = 1.d0
***
***********************************************************************
* The process is generated in the center of mass of the parton system
***********************************************************************
      p1(0) = x1pdf*ebeam1
      p2(0) = x2pdf*ebeam2
      if (p1(0).lt.m1.or.p2(0).lt.m2) then
         i_ok = 0
         return
      endif
      beta1 = 1.d0 !sqrt(1.d0 - m1**2/p1(0)**2)
      beta2 = 1.d0 !sqrt(1.d0 - m2**2/p2(0)**2)
      if (ipart1.eq.0.or.ipart2.eq.0) then
         beta1 = sqrt(1.d0 - m1**2/p1(0)**2)
         beta2 = sqrt(1.d0 - m2**2/p2(0)**2)
      endif

      p1(1) = 0.d0
      p1(2) = 0.d0
      p1(3) = beta1 * p1(0)

      p2(1) = 0.d0
      p2(2) = 0.d0
      p2(3) = -beta2 * p2(0)

      bis = (p1(3)+p2(3))/(p1(0)+p2(0))
      gis =  1.d0/sqrt(1.d0-bis**2)
      vboost(3) = bis

! "*cm" means CoM and up-type particle along +z!!!!!!!! 
      call boost(gis,vboost,p1,p1cm)
      call boost(gis,vboost,p2,p2cm)
*** 
* check if the event is a [O(alpha)] photon induced process
      ioalphind = 0
      if ((ipart1.eq.0.and.ipart2.ne.0).or.
     .    (ipart1.ne.0.and.ipart2.eq.0)) then
         ioalphind = 1

c SCOMMENTARE !!!!!
         if (ord.eq.'born') then
            i_ok = 0
            return
         endif
c SCOMMENTARE !!!!!

      endif
***
      if (ioalphind.eq.1) colfac = 1.d0

      imirror = 0
      if (boson.eq.'W') then
! I generate with up-type quark along +z
         if ((abs(ipart1).eq.2).or.(abs(ipart1).eq.4).or.ipart2.eq.0) 
     .        then
            continue
         else
            imirror = 1
         endif
      else ! if boson eq Z
         zIS = 'ddbar'
         if (mod(ipart1,2).eq.0) zIS = 'uubar'
         if (ioalphind.eq.0) then
            if (ipart1.ge.0) then
               continue
            else
               imirror = 1
            endif         
         else
            if (ipart1.eq.0) imirror = 1
         endif
      endif

*********************************************
**** QUANTITIES FOR up down SUBTRACTION CONES !!!
* and for cut efficiency...
cc      sqrtscalesub = p1cm(0) + p2cm(0)
c Maybe it is better to keep always mw for sqrtscalesub
      sqrtscalesub = spdf
      if (boson.eq.'W') then
         sqrtscalesub = mw
      else
         sqrtscalesub = mz
      endif
*************************************************
      if (imirror.eq.0) then
         ipart1mirr = ipart1
         ipart2mirr = ipart2
      else
!         p1cm and p2cm have to be mirrored and exchanged!
         ptmp(0) = p1cm(0)
         ptmp(1) = p1cm(1)
         ptmp(2) = p1cm(2)
         ptmp(3) = p1cm(3)
         p1cm(0) =  p2cm(0)
         p1cm(1) =  p2cm(1)
         p1cm(2) =  p2cm(2)
         p1cm(3) = -p2cm(3)
         p2cm(0) =  ptmp(0)
         p2cm(1) =  ptmp(1)
         p2cm(2) =  ptmp(2)
         p2cm(3) = -ptmp(3)
         ipart1mirr = ipart2
         ipart2mirr = ipart1
      endif

*** beginning event generation...
      einitial = p1cm(0)+p2cm(0)
      flux = 8.d0 * (einitial/2.d0)**2
ccc FIXED  in common    esoft = eps*einitial/2.d0
      if (i_pdf.eq.1) then
         esoftmax = 0.004d0
         if (einitial.gt.200.d0) esoftmax = 0.04d0
         if (einitial.gt.500.d0) esoftmax = 0.1d0
c!C!C!C!C!C!C
         esoftmax = eps*einitial/2.d0
c!C!C!C!C!C!C
      endif      
      
      esoft = esoftmax
      if (ioalphind.eq.0) then
         call multiplicityv2(eps,einitial,ng,wnphot)
      else
         ng = 1
         wnphot = 1.d0
      endif

      nsignal = 0
c      if (ng.lt.-1000) then
c         nsignal = ng
c         if (nsignal.eq.-1002) ng=0
c         if (nsignal.eq.-1003) ng=1
c      endif
      nphot = ng

      if (imirror.eq.0) then
         ama(1) = m1
         ama(2) = m2
      else
         ama(1) = m2
         ama(2) = m1
      endif

      if (boson.eq.'W') then
         ama(3) = mfs
         ama(4) = 0.d0
         ama3   = ama(3)
         charges(1) = 2.d0/3.d0
         charges(2) =-1.d0/3.d0
         charges(3) =-1.d0
         charges(4) = 0.d0
         amasses(1) = ama(1)
         amasses(2) = ama(2)
         amasses(3) = ama(3)
         amasses(4) = 0.d0
      else ! if Z
         ama(3) = mfs
         ama(4) = mfs
         ama3   = ama(3)

         carica = chq(abs(ipart1)) ! initialization
         if (abs(chq(ipart1)).gt.1.d-3) carica = chq(abs(ipart1)) 
         if (abs(chq(ipart2)).gt.1.d-3) carica = chq(abs(ipart2)) 

         charges(1) = carica
         charges(2) = charges(1)
         charges(3) =-1.d0
         charges(4) =-1.d0
         amasses(1) = ama(1)
         amasses(2) = ama(2)
         amasses(3) = ama(3)
         amasses(4) = ama(4)
      endif
****
      e1 = (einitial**2 + ama(1)**2 - ama(2)**2)/2.d0/einitial
      e2 = einitial - e1

      if (e1.lt.ama(1).or.e2.lt.ama(2)) then
         i_ok = 0
         return
      endif

      p1cm(0) = e1
      p1cm(1) = 0.d0
      p1cm(2) = 0.d0
      p1cm(3) = e1*sqrt(1.d0-ama(1)**2/e1**2)
      p2cm(0) =  e2
      p2cm(1) =  0.d0
      p2cm(2) =  0.d0
      p2cm(3) = -e2*sqrt(1.d0-ama(2)**2/e2**2)

      sdif = wnphot

      maaa = ama(1)
      mbbb = ama(2)

      if (ioalphind.eq.0) then
         call phasespace(p1cm,p2cm,p3cm,p4cm,qphcm,ng,
     .        ama3,ama(4),esoft,w,phsp,ie)
      else
         if (boson.eq.'W') then
            call phasespacephind(p1cm,p2cm,p3cm,p4cm,p5cm,qphcm,ama3,
     .           ama(4),ama(1),w,phsp,ie)
         else
            call phasespacephindz(p1cm,p2cm,p3cm,p4cm,p5cm,qphcm,ama3,
     .           ama(4),ama(1),w,phsp,ie)               
         endif
      endif
      
      if (ie.gt.1) ie = 1
      if (ie.gt.0) then
         i_ok = 0
         return
      endif
      
      wsampling  = w
      spaziofasi = phsp

      sdif = sdif * phsp * w
**********************************************
* Cuts are applied....
* all FS momenta are reported in the laboratory
      if (imirror.eq.1) then
         call mirror(p3cm,p3)
         call mirror(p4cm,p4)
         if (ioalphind.eq.1) call mirror(p5cm,p5)
      else
         do k = 0,3
            p3(k) = p3cm(k)
            p4(k) = p4cm(k)
            if (ioalphind.eq.1) p5(k) = p5cm(k)
         enddo
      endif
      vboost(3) = -bis
      call boost(gis,vboost,p3,p3)
      call boost(gis,vboost,p4,p4)
      if (ioalphind.eq.1) call boost(gis,vboost,p5,p5)
      if (ng.gt.0) then
         do kkk = 1,ng
            ptmp(0) = qphcm(kkk,0)
            ptmp(1) = qphcm(kkk,1)
            ptmp(2) = qphcm(kkk,2)
            ptmp(3) = qphcm(kkk,3)
            if (imirror.eq.1) call mirror(ptmp,ptmp)
            call boost(gis,vboost,ptmp,ptmp)
            qph(kkk,0) = ptmp(0) 
            qph(kkk,1) = ptmp(1) 
            qph(kkk,2) = ptmp(2)
            qph(kkk,3) = ptmp(3)
         enddo
      endif

      call cuts(p3,p4,qph,icut)

      if (icut.eq.1) then
         i_ok = 0
         return
      endif
**********************************************
      if ((boson.eq.'Z').and.(ischemealpha.ge.1)) then
         q2       = einitial*einitial
         dalpha   = 0.d0
         dalphaan = 0.d0
         if (ischemealpha.eq.2) then 
            dalpha   = 1.d0 - 1.d0/vpol(q2)
            dalphaan = 1.d0 - 1.d0/vpolan(q2)
         endif
         if (ifirst2.eq.0) then 
            drhoir1  = getdrhoir(1)
            drhoir12 = getdrhoir(2)
            rhofiu = getrhofi(mz*mz,1)
            rhofid = getrhofi(mz*mz,2)
            call getkappamuno(dkl,dku,dkd)
            if (ischemealpha.eq.1) then
ccc               print*,'Should not pass from here... event.f line ~398'
               dkl = 0.d0
               dku = 0.d0
               dkd = 0.d0
               drhoir1  = 0.d0
               drhoir12 = 0.d0
               rhofiu   = 1.d0
               rhofid   = 1.d0
            endif
         endif
         ifirst2 = 1

      endif

*** ng-photons real contribution
ccc      goto 111 ! phase space!!
      call real_photons(model,ng,einitial,p3cm,p4cm,qphcm,ie,
     >     imtx,emtx,iemap,emtxsub,ioalphind)
      
      if (iemap.gt.0) then
         i_ok = 0
         return
      endif
** soft+virtual contribution
      if (ie.eq.0) then
         sv = 1.d0
         if (ioalphind.eq.0) then
            if (nsignal.eq.0)         
     .           call svfactor(model,ng,einitial,p3cm,p4cm,eps,sv)
         endif
         sdif = sdif * sv
      endif
      svdebug = sv
****
      if (ng.eq.1.and.ioalphind.eq.0) then
         emtx = emtx + 1.d0 * emtxsub/sv
      endif

      emtx = emtx/in_conf_spin ! divided by initial spin conf. 
      emtx = emtx*convfac/flux ! divided by the flux
      emtx = emtx*Vud2*colfac
      sdif = sdif * emtx
      
      if (isnan(sdif)) then
         print*,'line 464 event.f'
         print*,'PID ',getpid()
         print*,boson,lepton,ord
         print*,sdif,sv,emtx,ng
         print*,'initial state',ipart1, ipart2
         print*,w,phsp
         print*,xxxx1,xxxx2
         print*,p1cm
         print*,p2cm
         print*,p3cm
         print*,p4cm        

         print*,'Resetting the event!!'
         sdif = 0.d0
      endif
****
      sdif = sdif * dpu
 111  continue ! phasespace!!

***SCOMMENTARE
      if ((isubtract.eq.1).and.ng.eq.0
     .     .and.i_pdf.eq.1.and.ie.eq.0
     .     .and.(ord.eq.'alpha'.or.ord.eq.'exp')) then

***SCOMMENTARE
***COMMENTARE
c      if ((isubtract.eq.1).and.ng.eq.0
c     .     .and.i_pdf.eq.1.and.ie.eq.0) then
***COMMENTARE

         sdborn = sdif/sv
        
         if (boson.eq.'Z') then
            sdborn = sdborn * bornmeclean/bornme
         endif

         call get_subtraction(x1pdf,x2pdf,bigq1,bigq2)
         q1 = pdf1(ipart1)
         q2 = pdf2(ipart2)
 
         if (ipart1.eq.0.and.ipart2.eq.0) then
! because the interpolation procedure returns always bigq < 0, and for
! gg it has to be > 0!
            bigq1 = -bigq1
            bigq2 = -bigq2
         endif

***SCOMMENTARE
         subtract = (bigq1/q1 + bigq2/q2 + subtractedoalsv)*sdborn
         sdif     = sdif - subtract
         
c         if (isnan(sdif)) then
c            print*,'event.f line 515'
c            print*,sdif,subtractedoalsv,sv
c         endif
***SCOMMENTARE
***COMMENTARE
c         subtract = (bigq1/q1 + bigq2/q2)*sdborn
c         sdif = - subtract *(-1.d0) !! NEEDED BECAUSE get_sub assumes
c                                   !! NEGATIVE VALUES FOR BIGQ!!!
c                                   !! WHICH IS NOT THE CASE FOR PHIND ALONE..
***COMMENTARE 
      endif

**********************************************
* The jacobian for the MC integration is calculated
      xjac = wwpwm
      xjac = xjac*wpdfpdfpdfjac          ! for pdfs...

*********************************************
      if (iscalechoice.gt.1) then     
***** PDF RUNNING SCALE *********************
         sdif = sdif/pdf1(ipart1)/pdf2(ipart2)
c     spdf = sqrt(shat)

         spdf = dot(p3cm+p4cm,p3cm+p4cm)
         if (iscalechoice.eq.2) then
            ! sqrt(M_leptons^2)
            spdf = sqrt(spdf)
         else
            ! sqrt(M_leptons^2 + pt_lepton_pair^2)
            spdf = sqrt(spdf+(p3cm(1)+p4cm(1))**2+(p3cm(2)+p4cm(2))**2)
         endif
         spdfcomm = spdf
         call get_pdf(x1pdf,x2pdf,spdf)
         sdif = sdif*pdf1(ipart1)*pdf2(ipart2)
*********************************************
*********************************************
      endif
      
      sdif = sdif*xjac
      sumsez(ng+1)  = sumsez(ng+1) + sdif
      sumsez2(ng+1) = sumsez2(ng+1)+ sdif*sdif
      sezij(ipart1,ipart2) = sezij(ipart1,ipart2) + sdif

      sezij2(ipart1,ipart2) = sezij2(ipart1,ipart2) + sdif*sdif

      icounterpartial(ng+1) = icounterpartial(ng+1) + 1
      if (ng.gt.nphotmax) nphotmax = ng

**************************************************
** assigning ipart3,ipart4,ipart5 MC id according to 
** PDG MC scheme, simplyfing event storage to be processed with
** Shower MCs. They correspond to momenta p3,p4,qph(1,:)
      imakemccodesavailable = 1
      if (storing.eq.'y'.or.imakemccodesavailable.eq.1) then
         ipart5 = -222222222
         if (ioalphind.eq.0.and.ng.gt.0) ipart5 = 22
         if (ioalphind.eq.1) then
            if (boson.eq.'Z') then
               ipart5 = ipart2
               if (ipart1.ne.0) ipart5 = ipart1
            else ! if W
               !ckm-ology !!
               call wraprng(csi,1)
** d s b -  1 3 5
** u c t -  2 4 6
               ibs1 = abs(ipart1)
               ibs2 = abs(ipart2)
               isg1 = ipart1/max(0.1,1.*ibs1)
               isg2 = ipart2/max(0.1,1.*ibs2)
               if (ibs1.eq.1.or.ibs1.eq.3.or.ibs1.eq.5) then
                  p_1 = myckm(ipart1,2)**2
                  p_2 = myckm(ipart1,4)**2
                  p_3 = 0.d0
                  id1 = 2*isg1
                  id2 = 4*isg1
                  id3 = 0*isg1
               endif
               if (ibs1.eq.2.or.ibs1.eq.4) then
                  p_1 = myckm(ipart1,1)**2
                  p_2 = myckm(ipart1,3)**2
                  p_3 = myckm(ipart1,5)**2
                  id1 = 1*isg1
                  id2 = 3*isg1
                  id3 = 5*isg1
               endif
               if (ibs2.eq.1.or.ibs2.eq.3.or.ibs2.eq.5) then
                  p_1 = myckm(ipart2,2)**2
                  p_2 = myckm(ipart2,4)**2
                  p_3 = 0.d0
                  id1 = 2*isg2
                  id2 = 4*isg2
                  id3 = 0*isg2
               endif
               if (ibs2.eq.2.or.ibs2.eq.4) then
                  p_1 = myckm(ipart2,1)**2
                  p_2 = myckm(ipart2,3)**2
                  p_3 = myckm(ipart2,5)**2
                  id1 = 1*isg2
                  id2 = 3*isg2
                  id3 = 5*isg2
               endif
               s_pi = p_1+p_2+p_3
               p_1 = p_1/s_pi
               p_2 = p_2/s_pi
               p_3 = p_3/s_pi
               if (csi(1).le.p_1) then
                  ipart5 = id1
               elseif (csi(1).le.(p_1+p_2)) then
                  ipart5 = id2
               else
                  ipart5 = id3
               endif
            endif
         endif
         if (boson.eq.'W') then
            if (lepton.eq.1) then
               if (chfs.gt.0.d0) then
                  ipart3 = -11
                  ipart4 =  12
               else
                  ipart3 =  11
                  ipart4 = -12
               endif
            endif
            if (lepton.eq.2) then
               if (chfs.gt.0.d0) then
                  ipart3 = -13
                  ipart4 = 14
               else
                  ipart3 = 13
                  ipart4 = -14
               endif
            endif
            if (lepton.eq.3) then
               if (chfs.gt.0.d0) then
                  ipart3 = -15
                  ipart4 = 16
               else
                  ipart3 = 15
                  ipart4 = -16
               endif
            endif
         else ! Z
            if (lepton.eq.1) then
               ipart3 =  11
               ipart4 = -11
            elseif (lepton.eq.2) then
               ipart3 =  13
               ipart4 = -13
            else
               ipart3 =  15
               ipart4 = -15
            endif
         endif
         ip1mc = ipart1
         ip2mc = ipart2
         ip3mc = ipart3
         ip4mc = ipart4
         ip5mc = ipart5
         if (ip1mc.eq.0) ip1mc = 22
         if (ip2mc.eq.0) ip2mc = 22
      endif
**************************************************
      return
      end
********************************************************
      subroutine generate_event_OLDPS(sdif,i_ok,p3,p4,qph)
      include 'shared.inc'
      common/n_phot/nphot

      real*4 csi(1),rnd(1)
      dimension x(40),xvect(10),vboost(0:3)
      dimension qphtmp(40,0:3)
      common/firstgenevent/ifirst
      data ifirst/1/
      do i=0,3
         vboost(i) = 0d0
      enddo

      nphot = 0
      
      if (ifirst.eq.0) then
         do i = 1,20
            sumsez(i)   = 0.d0
            sumsez2(i)  = 0.d0
            sezioni(i)  = 0.d0
            sezioni2(i) = 0.d0
            icounterpartial(i) = 0
         enddo
      endif

      i_ok = 1
      s = 4.d0*ebeam1*ebeam1
      call get_x1x2pdf(s,pdfjac)
      
      spdf = sqrt(4.d0*x1pdf*x2pdf*ebeam1*ebeam2) ! sqrt of the pdf Q^2 scale

      spdf = mw * rescalepdfscale
      if (boson.eq.'Z') spdf = mz * rescalepdfscale

      spdfcomm = spdf
   
      if (spdf.lt.qpdfmin) then
         i_ok = 0
         return 
      endif

************* W+ W- production at the same time!!
      wwpwm = 1.d0
      if (abs(chfs).gt.0.d0.and.ichfs.eq.2) then
         chfs = 1.d0
         call wraprng(csi,1)
         if (csi(1).le.0.5) chfs = -1.d0
         wwpwm = 2.d0
      endif
**************************************
      call get_is(x1pdf,x2pdf,spdf,wpdfpdf)
      if (ipart1.eq.(ifl+1)) then
         i_ok = 0
         return
      endif
***********************************************************************
* The process is generated in the center of mass of the parton system
*
* The cos_theta and phi integration variable are extracted
* ccm is sampled for W according (1-c)**2 or (1+c)**2
      call wraprng(rnd,1)
      if (boson.eq.'W') then
         if ((abs(ipart1).eq.2).or.(abs(ipart1).eq.4)) then   
            ccm = 1.d0-2.d0*(1.d0-rnd(1))**(1.d0/3.d0)
            xjacccm = 8.d0/3.d0/(1.d0-ccm)/(1.d0-ccm)         
         else
            ccm = 2.d0*(1.d0-rnd(1))**(1.d0/3.d0)-1.d0
            xjacccm = 8.d0/3.d0/(1.d0+ccm)/(1.d0+ccm)
         endif 
      else  ! if (boson.eq.'Z')
c      ccm = 2.d0 * rnd(1) - 1.d0
c      xjacccm = 2.d0
         ccm = ccmgen(wei)
         xjacccm = wei
      endif
      call wraprng(rnd,1)
      phi  = 2.d0*pi*rnd(1)

      bis = (x1pdf*ebeam1-x2pdf*ebeam2)/(x1pdf*ebeam1+x2pdf*ebeam2)
                                      ! beta (speed) of the center of mass
      gis =  1.d0/sqrt(1.d0-bis**2)   ! and gamma....

      E1l = x1pdf*ebeam1
      E2l = x2pdf*ebeam2
* the energy in the parton CoM is calculated. Of course, E1 and E2 are equal!
      E1 = gis*E1l*(1.d0-bis)
      E2 = gis*E2l*(1.d0+bis)

! 18 / 8 / 2003
!!! I recalculate E1 and E2 for massive quarks. I use the convention that
!!! xpdf are fraction of proton p_z carried by quarks, as agreed with Wieslaw.
!!! However, in my opinion, the ambiguity remains....

      if (i_scheme.eq.1) then 
!! I use the "p_z scheme"....
         pz1l =  x1pdf*ebeam1
         pz2l = -x2pdf*ebeam2   ! it is along the negative z axis
      
         E1l = sqrt(pz1l**2 + m1**2)
         E2l = sqrt(pz2l**2 + m2**2)
      endif
!! or the "E scheme"....
      if (i_scheme.eq.2) then
         E1l = x1pdf*ebeam1
         E2l = x2pdf*ebeam2
      
         if (E1l.lt.m1.or.E2l.lt.m2) then
            i_ok = 0
            return
         endif
         pz1l = sqrt(E1l**2 - m1**2)
         pz2l = -sqrt(E2l**2 - m2**2)
      endif
!!
!! or the "light cone scheme"....
      if (i_scheme.eq.3) then
         sqr2 = sqrt(2.d0)
         
         plcp1 = 2.d0 * ebeam1 / sqr2
         plcp2 = 2.d0 * ebeam2 / sqr2
         
         qlc1 =  x1pdf * plcp1
         qlc2 =  x2pdf * plcp2

         E1l  = 2.d0*qlc1**2+m1**2
         E1l  = E1l/2.d0/sqr2/qlc1
         
         E2l  = 2.d0*qlc2**2+m2**2
         E2l  = E2l/2.d0/sqr2/qlc2
         
         if (E1l.lt.m1.or.E2l.lt.m2) then
            i_ok = 0
            return
         endif
         
         pz1l = sqrt(E1l**2-m1**2)
         pz2l = -sqrt(E2l**2-m2**2)
      endif
!!
      if (i_scheme.eq.4) then
         E1l = x1pdf*ebeam1
         E2l = x2pdf*ebeam2
         
         pz1l =  E1l
         pz2l = -E2l         
      endif

      E1  = gis*(E1l - pz1l*bis)
      E2  = gis*(E2l - pz2l*bis)

      E1cm = E1 ! to be passed in a common, since E1 and E2 conflict with a  
      E2cm = E2 ! soubroutine argument....

      pz1  = gis*(pz1l - E1l*bis)
      pz2  = gis*(pz2l - E2l*bis)
*
* the QED SF scale is set...
      q2 = 4.d0*(E1)**2
* .... and the QED PS is called for the final state particles
      call qed_stuff(E1,E2,q2,ccm,phi,p3,p4,qph,x1,x2,x3,x4,ierror)
      if (ierror.gt.0) then
         i_ok = 0
         return
      endif

********************************
********************************
      ! compact qph at the beginning of the array...
      do k = 0,3
         do j = 1,40
            qphtmp(j,k) = 0.d0
         enddo
      enddo
      i = 1
      do k = 1,40
         if (qph(k,0).gt.0.d0) then
            nphot = nphot + 1
            qphtmp(i,0) = qph(k,0)
            qphtmp(i,1) = qph(k,1)
            qphtmp(i,2) = qph(k,2)
            qphtmp(i,3) = qph(k,3)
            qph(k,0) = 0.d0
            qph(k,1) = 0.d0
            qph(k,2) = 0.d0
            qph(k,3) = 0.d0
            i = i + 1
         endif
      enddo
      do k = 0,3
         do j = 1,40
            qph(j,k) = qphtmp(j,k)
         enddo
      enddo
********************************
********************************

* the 'shat' of the process is calculated
      shat_old = 4.d0*x1*x2*E1*E2

! 19 / 8 / 2003 ! I calculate massive shat
      shat = m1**2 + m2**2 + 2.d0 * (E1*E2 - pz1 * pz2)
      if (i_scheme.eq.4) shat = 2.d0 * (E1*E2 - pz1 * pz2)
      ! which is the same as:
      !shataaa = m1**2 + m2**2 + 2.d0 * (E1l*E2l - pz1l * pz2l)
**********************************************
* Cuts are applied.....
* all FS momenta are reported in the laboratory
      vboost(3) = -bis
      call boost(gis,vboost,p3,p3)
      call boost(gis,vboost,p4,p4)      
      ng = 0
      do kkk = 1,40
         if (qph(kkk,0).gt.0.d0) then
            ng = ng + 1
            ptmp(0) = qph(kkk,0)
            ptmp(1) = qph(kkk,1)
            ptmp(2) = qph(kkk,2)
            ptmp(3) = qph(kkk,3)
            call boost(gis,vboost,ptmp,ptmp)
            qph(kkk,0) = ptmp(0) 
            qph(kkk,1) = ptmp(1) 
            qph(kkk,2) = ptmp(2)
            qph(kkk,3) = ptmp(3)         
         endif
      enddo
      call cuts(p3,p4,qph,icut)
      if (icut.eq.1) then
         i_ok = 0
         return
      endif
********************************************
      i1 = ipart1
      i2 = ipart2
      sdif = shard(shat,ccm,i1,i2)
* The jacobian for the MC integration is calculated
!      xjac = (xmax-xmin)*(xmax-xmin) ! this is for x1-2pdf integration
      xjac = pdfjac       * wwpwm
      xjac = xjac*2.d0*pi*xjacccm  ! this is for ccm and phi integration
      xjac = xjac*wpdfpdf          ! for pdfs...
      sdif = sdif*xjac
      ifirst = 0


*********************************************
*********************************************
      if (iscalechoice.gt.1) then     
***** PDF RUNNING SCALE *********************
         sdif = sdif/pdf1(ipart1)/pdf2(ipart2)
c     spdf = sqrt(shat)

         spdf = dot(p3+p4,p3+p4)
         if (iscalechoice.eq.2) then
            ! sqrt(M_leptons^2)
            spdf = sqrt(spdf)
         else
            ! sqrt(M_leptons^2 + pt_lepton_pair^2)
            spdf = sqrt(spdf+(p3(1)+p4(1))**2+(p3(2)+p4(2))**2)
         endif
         spdfcomm = spdf
         call get_pdf(x1pdf,x2pdf,spdf)
         sdif = sdif*pdf1(ipart1)*pdf2(ipart2)
*********************************************
*********************************************
      endif

      
      sezij(ipart1,ipart2) = sezij(ipart1,ipart2) + sdif
**************************************************
** assigning ipart3,ipart4,ipart5 MC id according to 
** PDG MC scheme, simplyfing event storage to be processed with
** Shower MCs. They correspond to momenta p3,p4,qph(1,:)!
      imakemccodesavailable = 1
      if (storing.eq.'y'.or.imakemccodesavailable.eq.1) then
         ioalphind = 0
         ipart5 = -222222222
         if (ioalphind.eq.0.and.ng.gt.0) ipart5 = 22
         if (ioalphind.eq.1) then
            if (boson.eq.'Z') then
               ipart5 = ipart2
               if (ipart1.ne.0) ipart5 = ipart1
            else ! if W
               !ckm-ology !!
               call wraprng(csi,1)
** d s b -  1 3 5
** u c t -  2 4 6
               ibs1 = abs(ipart1)
               ibs2 = abs(ipart2)
               isg1 = ipart1/max(0.1,1.*ibs1)
               isg2 = ipart2/max(0.1,1.*ibs2)
               if (ibs1.eq.1.or.ibs1.eq.3.or.ibs1.eq.5) then
                  p_1 = myckm(ipart1,2)**2
                  p_2 = myckm(ipart1,4)**2
                  p_3 = 0.d0
                  id1 = 2*isg1
                  id2 = 4*isg1
                  id3 = 0*isg1
               endif
               if (ibs1.eq.2.or.ibs1.eq.4) then
                  p_1 = myckm(ipart1,1)**2
                  p_2 = myckm(ipart1,3)**2
                  p_3 = myckm(ipart1,5)**2
                  id1 = 1*isg1
                  id2 = 3*isg1
                  id3 = 5*isg1
               endif
               if (ibs2.eq.1.or.ibs2.eq.3.or.ibs2.eq.5) then
                  p_1 = myckm(ipart2,2)**2
                  p_2 = myckm(ipart2,4)**2
                  p_3 = 0.d0
                  id1 = 2*isg2
                  id2 = 4*isg2
                  id3 = 0*isg2
               endif
               if (ibs2.eq.2.or.ibs2.eq.4) then
                  p_1 = myckm(ipart2,1)**2
                  p_2 = myckm(ipart2,3)**2
                  p_3 = myckm(ipart2,5)**2
                  id1 = 1*isg2
                  id2 = 3*isg2
                  id3 = 5*isg2
               endif
               s_pi = p_1+p_2+p_3
               p_1 = p_1/s_pi
               p_2 = p_2/s_pi
               p_3 = p_3/s_pi
               if (csi(1).le.p_1) then
                  ipart5 = id1
               elseif (csi(1).le.(p_1+p_2)) then
                  ipart5 = id2
               else
                  ipart5 = id3
               endif
            endif
         endif
         if (boson.eq.'W') then
            if (lepton.eq.1) then
               if (chfs.gt.0.d0) then
                  ipart3 = -11
                  ipart4 = 12
               else
                  ipart3 = 11
                  ipart4 = -12
               endif
            endif
            if (lepton.eq.2) then
               if (chfs.gt.0.d0) then
                  ipart3 = -13
                  ipart4 = 14
               else
                  ipart3 = 13
                  ipart4 = -14
               endif
            endif
            if (lepton.eq.3) then
               if (chfs.gt.0.d0) then
                  ipart3 = -15
                  ipart4 = 16
               else
                  ipart3 = 15
                  ipart4 = -16
               endif
            endif
         else ! Z
            if (lepton.eq.1) then
               ipart3 =  11
               ipart4 = -11
            elseif (lepton.eq.2) then
               ipart3 =  13
               ipart4 = -13
            else
               ipart3 =  15
               ipart4 = -15
            endif
         endif
         ip1mc = ipart1
         ip2mc = ipart2
         ip3mc = ipart3
         ip4mc = ipart4
         ip5mc = ipart5
         if (ip1mc.eq.0) ip1mc = 22
         if (ip2mc.eq.0) ip2mc = 22
      endif
**************************************************
      return
      end
********************************************************************
      function shard(s,c,ip,jp)
      include 'shared.inc'
      double precision mz2
      E1 = E1cm
      E2 = E2cm
      if (boson.eq.'W') then
*************** SHARD for W production
* c is the cosine of the angle between the up-type quark and the lepton.
* In any case this function is symmetric under c --> -c
         Vud2   = 1.d0 !! now it is in the PDFs !!!
         s4tw   = s2tw**2
         
         gw_run = gw*s/mw**2
         gw_run = gw            ! fixed W width !!!
c         call messages(7)

         if ((ip.eq.2).or.(ip.eq.4))   then 
            coseno =  c
            qma    =  m1
            pquark =  pz1
            equark =  E1
         endif
         if ((ip.eq.-1).or.(ip.eq.-3)) then 
            coseno = -c
            qma    =  m2
            pquark =  -pz2
            equark =  E2
         endif
         if ((ip.eq.-2).or.(ip.eq.-4)) then
            coseno =  c
            qma    =  m1
            pquark =  pz1
            equark =  E1
         endif
         if ((ip.eq.1).or.(ip.eq.3)) then
            coseno = -c 
            qma    =  m2
            pquark =  -pz2
            equark =  E2
         endif
         if (ip.eq.-5) then 
            coseno = -c
            qma    =  m2
            pquark =  -pz2
            equark =  E2
         endif
         if (ip.eq.5) then
            coseno = -c
            qma    =  m2
            pquark =  -pz2
            equark =  E2
         endif
         if (i_scheme.eq.4) qma = 0.d0
!!!!!!!!!!         if ((ip.eq.-2).or.(ip.eq.-4)) coseno = -c
!!!!!!!!!!         if ((ip.eq.1).or.(ip.eq.3))   coseno =  c
!! I take massless final state lepton....
         plepton = (E1 + E2)/2.d0
         elepton = plepton

         u_man = qma**2-2.d0*(equark*elepton-coseno*pquark*plepton)
         u2    = u_man**2
         
         u_man_old  = -s/2.d0*(1.d0-coseno)
         u2_old     = u_man_old**2

         den    = (s-mw**2)**2+mw**2*gw_run**2         
         alp = alpha_gf
         if (ischemealpha.eq.0) alp = alpha

         shard  = alp**2*Vud2/s4tw/48.d0 ! /192.d0
         shard  = shard/s*u2
         shard  = shard/den
         shard  = shard * convfac
      else ! if (boson.eq.'Z') then
*************** SHARD for Z production
         modpart = abs(ip)
         kpart   = modpart
         qq      = chq(kpart)

         spine = -0.5d0  ! weak electron isospin
         qe    = -1.d0

         gvq = itre(kpart) - 2.d0*qq*s2tw 
         gaq = itre(kpart)

         gve = spine - 2.d0*qe*s2tw
         gae = spine

         gz_run = gz*s/mz**2
         gz_run = gz      ! fixed Z width !!!
         mz2    =  mz**2
         den    = (s-mz2)**2+mz2*gz_run**2
         if (ip.gt.0) coseno =  c
         if (ip.lt.0) coseno = -c
         r2 = sqrt(2.d0)

!         goto 333
c$$$c         if (ialphachoice.eq.0) alp = alpha
c$$$c         if (ialphachoice.eq.1) alp = alpha_gf
c$$$c         if (ialphachoice.eq.2) alp = alpha*vpol(s)
c$$$c         if (ialphachoice.ge.3) call messages(3)
c$$$         alp = alpha
c$$$         sgffa  = s*gf/alp
c$$$         sgffa2 = sgffa**2
c$$$         qq2 = qq**2
c$$$         fs = qq2 + qe*qq*(gve*gvq/pi/r2 * mz2*(s-mz2)/den)*sgffa
c$$$         fs = fs + (gve**2+gae**2)*(gvq**2+gaq**2)/8/pi/pi *
c$$$     >        mz2**2/den*sgffa2
c$$$         gs = qe*qq*r2*gae*gaq/pi*mz2*(s-mz2)/den*sgffa
c$$$         gs = gs + gve*gae*gvq*gaq/pi/pi*mz2**2/den*sgffa2
c$$$         shard = alp**2/4.d0/s*(fs*(1.d0+coseno**2) + gs*coseno)         
c$$$         shard = shard * convfac * 1.d0 / 3.d0 !! color factor....
c$$$         return
*----------------------------------------------------
! more clear to implement the EW scheme....        
 333     continue

         if (ischemealpha.eq.0) then 
            alp = alpha
            e2 = 4.d0*pi*alp 
            e4 = e2*e2
            gfp = r2*e2/s2tw/8.d0/mw/mw ! ALPHA(0) SCHEME!!!!!!!
         else
            alp = alpha*vpol(s)
            e2 = 4.d0*pi*alp 
            e4 = e2*e2
            gfp = gf
            e2p = 4.d0*pi*alpha_gf
            gfp = r2*e2p/s2tw/8.d0/mw/mw
         endif

         modpart = abs(ip)
         qq      = chq(modpart)
         qe      = -1.d0
         spine   = -0.5d0

         gvq = itre(modpart) - 2.d0*qq*s2tw 
         gaq = itre(modpart)

         gve = spine - qe * 2.d0*s2tw
         gae = spine

         fatt  = gfp*mz2*r2             ! prop. to (g/cos_thw)^2 !!
         rechi = fatt * s * (s-mz2)/den ! "  Z diagram    "
         chim2 = fatt**2 * s**2 / den   ! " |Z diagram|^2 "
         a0 = qq**2*qe**2*e4 *1.d0      !   |photon|^2 
         a0 = a0 +
     >        2.d0*gvq*gve*qq*qe*e2*rechi *1.d0     *1.d0 +! photon*Z
     >        (gvq**2+gaq**2)*(gve**2+gae**2)*chim2 *1.d0  ! |Z|^2
         a1 = 4.d0*qq*qe*gaq*gae*e2*rechi *1.d0     *1.d0 +! photon*Z
     >        8.d0*gvq*gaq*gve*gae*chim2            *1.d0  ! |Z|^2

         sss   = a0*(1.d0 + coseno**2) + a1 * coseno
         sss   = sss/64.d0/pi/pi/s * convfac /3.d0 
         shard = sss
*----------------------------------------------------
      endif
      end
