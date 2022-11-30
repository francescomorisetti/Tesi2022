**********************************************
      subroutine initrng(iseed)
! it can be modified by the user
      integer iseed,lux,k1,k2
      parameter (lux = 4,k1 = 0,k2 = 0) ! for RANLUX
      call rluxgo(lux,iseed,k1,k2)
      return
      end
*
      subroutine wraprng(csi,n)
! it can be modified by the user
      real*4 csi(n)
      call ranlux(csi,n)
      return
      end
**********************************************
      subroutine fillfeynartscommon
* this routine fills the common needed for FeynArts virtual libraries
      character*1 boson
      common/vectorboson/boson
      double precision EL, Alfa, Alfa2, GF, GS, Alfas, Alfas2
      double precision MW, MW2, MZ, MZ2, GW, GZ
      double precision SW, SW2, CW, CW2
      double precision MH, MH2, MG0, MG02, MGp, MGp2
      double precision ME, ME2, MM, MM2, ML, ML2, MLE(3), MLE2(3)
      double precision MU, MU2, MC, MC2, MT, MT2, MQU(3), MQU2(3)
      double precision MD, MD2, MS, MS2, MB, MB2, MQD(3), MQD2(3)
      double precision m1h,m2h,ch1h,ch2h,chfsh,mfsh
      integer leptonh,imirrorh
      double precision meh,mmuh,mwh,mzh,muh,mdh,mch,msh,mtauh,mbh,mth,
     .     mqh,chqh,itreh,mhh,gwh,gzh,
     .     alphah,alpha_gfh,s2twh,gfh,pih,convfach
      common/feynartshorace/Alfa, Alfa2, GF,
     .     MW, MW2, MZ, MZ2, GW, GZ, SW, SW2, CW, CW2,
     .     MH, MH2, ME, ME2, MM, MM2, ML, ML2, 
     .     MU, MU2, MC, MC2, MT, MT2, MD, MD2, MS, MS2, MB, MB2
      common/processhor/m1h,m2h,ch1h,ch2h,chfsh,mfsh,leptonh,imirrorh
      common/masseshor/meh,mmuh,mwh,mzh,muh,mdh,mch,msh,mtauh,mbh,mth,
     .     mqh(-5:5),chqh(-5:5),itreh(-5:5)
      common/higgsmass/mhh
      common/t_widths/gwh,gzh
      common/const/alphah,alpha_gfh,s2twh,gfh,pih,convfach
      double precision etamaxh,ptminh,cminh,cmaxh,ptminmissh,
     .     tmminh,mass2h,dalpha,dalphaan,drho,drhoir1,drhoir12,rhofiu,
     .     rhofid,dkl,dku,dkd
      integer ismearh,irech,icutmth,ischemealphah,ischeme
      common/exp_cuts/etamaxh,ptminh,cminh,cmaxh,ptminmissh,
     .     tmminh,mass2h,ismearh,irech,icutmth,ischemealphah
      common/feyhorflags/dalpha,dalphaan,drho,drhoir1,drhoir12,rhofiu,
     .     rhofid,dkl,dku,dkd,ischeme

      alfa  = alphah
      alfa2 = alfa*alfa
      gf    = gfh
      mw    = mwh
      mw2   = mw*mw
      mz    = mzh
      mz2   = mz*mz
      gw    = gwh
      gz    = gzh
      sw2   = s2twh
      sw    = sqrt(sw2)
      cw    = mw/mz
      cw2   = cw*cw
      mh    = mhh
      mh2   = mh*mh

      ischeme = ischemealphah

****      
      if (boson.eq.'W') then
! FA library was produced originally for e+ final state in the W case!!!
         if (leptonh.eq.1) then 
            me  = meh
            me2 = me*me
            mm  = mmuh 
            mm2 = mm*mm
            ml  = mtauh
            ml2 = ml*ml
         elseif (leptonh.eq.2) then
            mm  = meh
            mm2 = mm*mm
            me  = mmuh 
            me2 = me*me
            ml  = mtauh
            ml2 = ml*ml
         else
            ml  = meh
            ml2 = ml*ml
            me  = mtauh 
            me2 = me*me
            mm  = mmuh
            mm2 = mm*mm
         endif
      else ! if boson = 'Z'
! FA library was produced originally for mu final state in the Z case!!!
         if (leptonh.eq.1) then 
            mm  = meh
            mm2 = mm*mm
            me  = mmuh 
            me2 = me*me
            ml  = mtauh
            ml2 = ml*ml
         elseif (leptonh.eq.2) then
            me  = meh
            me2 = me*me
            mm  = mmuh 
            mm2 = mm*mm
            ml  = mtauh
            ml2 = ml*ml
         else
            mm  = mtauh
            mm2 = mm*mm
            me  = meh 
            me2 = me*me
            ml  = mmuh
            ml2 = ml*ml
         endif
      endif
****
      mu  = muh
      mu2 = mu*mu
      md  = mdh
      md2 = md*md
      mc  = mch
      mc2 = mc*mc
      ms  = msh
      ms2 = ms*ms
      mb  = mbh
      mb2 = mb*mb
      mt  = mth
      mt2 = mt*mt
      return
      end
************************************************************************
      subroutine print_output(iev,xsect,var,
     <     wsect,wvar2,eff,i_over,nneg,wn,wnv2,wo,wvo2,ncut)
      include 'shared.inc'
      integer*8 ncut,nwei,nunwneg
      common/nevents/nwei,nunwneg
      character*30 whichdeltar
      common/which_dr/whichdeltar
      common/deltaerre/deltar
      integer*8 incalls
      real*4 etres,etres12(2)
      character*6 ord
      common/qedORDER/ord
      common/tetto/tetto
      character*20 machinename
      character*20 who
      integer*8 istarttime
      common/starttime/istarttime

      open(10,file=outfile,status='unknown')
      rewind(10)
      rncalls = iev/eff
      callsn  = iev/eff
      cuteff  = 1.d0 - (1.d0*nwei)/(callsn)
      write(10,*)' '
      if (idumpweighted.eq.1)
     .     write(10,*)'Requested number of weighed events: ',nmax
      if (idumpweighted.eq.0)
     .     write(10,*)'Requested number of unweighed events: ',nmax
      write(10,*)'Unweighted event n. ',iev
      write(10,*)'Weighted event n. ',rncalls
      write(10,*)'Weighted events after cuts ',nwei
      write(10,*)' '
      if (i_pdf.eq.1) then
         write(10,'(1x,A,I5,2x,I5)')
     >        'Initial state particles: ',hadr1,hadr2
         write(10,'(1x,A,A)')
     >        'PDFs: ',commentonpdfs
      write(10,*)' '
      endif
      if (i_pdf.ne.1) write(10,'(1x,A,I4,I4)')
     >     'Initial state particles: ',iquark1,iquark2
      write(10,'(1x,A,f10.2,A)')'C. of m. energy =',2.d0*ebeam,' GeV'
      if (lepton.eq.1.and.ichfs.eq.-1) 
     .     write(10,*)'W- production, decaying into electrons'
      if (lepton.eq.1.and.ichfs.eq.1) 
     .     write(10,*)'W+ production, decaying into electrons'
      if (lepton.eq.1.and.ichfs.eq.2) 
     .     write(10,*)'W+ & W- production, decaying into electrons'
      if (lepton.eq.1.and.ichfs.eq.0) 
     .     write(10,*)'Z/photon production, decaying into electrons'
      if (lepton.eq.2.and.ichfs.eq.-1)
     .     write(10,*)'W- production, decaying into muons'
      if (lepton.eq.2.and.ichfs.eq.1) 
     .     write(10,*)'W+ production, decaying into muons'
      if (lepton.eq.2.and.ichfs.eq.2) 
     .     write(10,*)'W+ & W- production, decaying into muons'
      if (lepton.eq.2.and.ichfs.eq.0) 
     .     write(10,*)'Z/photon production, decaying into muons'
      if (lepton.eq.3.and.ichfs.eq.-1)
     .     write(10,*)'W- production, decaying into taus'
      if (lepton.eq.3.and.ichfs.eq.1) 
     .     write(10,*)'W+ production, decaying into taus'
      if (lepton.eq.3.and.ichfs.eq.2) 
     .     write(10,*)'W+ & W- production, decaying into taus'
      if (lepton.eq.3.and.ichfs.eq.0) 
     .     write(10,*)'Z/photon production, decaying into taus'


      write(10,*)' '
      write(10,*)'W mass      = ',mw,' GeV'
      write(10,*)'W width     = ',gw,' GeV'
      write(10,*)'Z mass      = ',mz,' GeV'
      write(10,*)'Z width     = ',gz,' GeV'
      write(10,*)'sin^2(th_w) = ',s2tw
      write(10,*)'delta_r     = ',deltar,whichdeltar
      write(10,*)' '
      if (iewk.eq.0) write(10,*)'QED RC: ',ord,' (OLD QED PS)'
      if (iewk.eq.1) write(10,*)'EW RC: ',ord,' (MATCHED QED PS)'
c      write(10,*)'Smearing and recombining flags: ',ismear,irec
      isca = ischemealpha
      if (isca.eq.2) isca = 1
      write(10,*)' EW input scheme: ',isca
      write(10,*)' '
      write(10,'(1x,A,f16.7,A,f16.7,A)')'hit or miss cross section  = ',
     >     xsect,' +- ',var,' (pb)'
                            ! 1x is a space!!!!!
      write(10,'(1x,A,f16.7,A,f16.7,A)')
     >     'h. or m. out-of-range bias = ',
     >     wo,' +-',sqrt(abs(wvo2)), ' (pb)'
      arg = dabs(wvar2)
      sw2 = sqrt(arg)
      write(10,'(1x,A,f16.7,A,f16.7,A)')'weighted cross section     = ',
     >     wsect,' +- ',sw2,' (pb)'
      write(10,*)' '
      write(10,'(1x,A,f9.5,A,I6,A)')'h. or m. efficiency = ',100*eff,
     >     ' %    (out of range = ',i_over,')'
      tobeset = ratio_max * tune_hm
c      write(10,*)'ratio max =',ratio_max,' (upper limit =',tune_hm,
c     >', to be set to',tobeset,')'
      write(10,*)'(upper limit =',tune_hm,', to be set to',tobeset,')'

      if (storing.eq.'y') then
ccc   CONFIGURATION FILE FOR HORACE LHA INTERFACE
         xmaxup = tetto/tune_hm*tobeset
         call writeconffile(wsect,sw2,xmaxup)
      endif
ccc

      write(10,*)' '
      if (storing.eq.'y')  then 
         write(10,*)'Events stored in file: ',storfile
      else
         write(10,*)'Events not stored'
      endif

      somma = 0.d0
      if (iewk.eq.1) then
      write(10,*)'Partial cross sections:'
      do i = 0,nphotmax
         ii = icounterpartial(i+1)
         seztot(i+1) = sumsez(i+1)/callsn
         sezerr(i+1) = abs(sumsez2(i+1)/callsn - seztot(i+1)**2)
         sezerr(i+1) = sqrt(sezerr(i+1)/callsn)
         write(10,*)' ',i,' photons:',seztot(i+1),' +-',sezerr(i+1),'
     .        (npoints: ',ii,')'
         somma = somma + seztot(i+1)
      enddo
      write(10,*)' '
      endif
      write(10,*)'Negative weights statistics:'
      write(10,*)' n. of weighted points   =',nneg
      write(10,*)' n. of unweighted events =',nunwneg
      write(10,*)' cross section ',wn,' +-',sqrt(abs(wnv2))
      write(10,*)' '
c      write(10,*)'cut points: ',cuteff*100.d0,' %'
      write(10,*)'partial q-q[''][bar] contributions'
      sijsum = 0.d0
      do i = -ifl,ifl
         do j = -ifl,ifl
            sij = sezij(i,j)/callsn
            if (abs(sij).gt.0.d0) then

               sije = abs(sezij2(i,j)/callsn - sij**2)
               sije = sqrt(sije/callsn)
               
               
               sijsum = sijsum + sij
               ii = i
               jj = j
               if (ii.eq.0) ii = 22
               if (jj.eq.0) jj = 22
             write(10,*)ii,jj,' ',sij,'+-',sije,'(',sije/sij*100d0,'%)'
            endif
         enddo
      enddo
      write(10,*)'Sum: ',sijsum
***
** hostname & user ID stuff!!
      ihn = hostnm(machinename)
      ipid = getpid()
      iuid = getuid()
      ilm = lnblnk(machinename)
***
      write(10,*)'....'
      write(10,*)'Started by user ID',iuid,' on ',machinename(1:ilm),
     .', process ID',ipid
      write(10,*)'....'
      call Etime(etres12,etres)
      write(10,*)'Running time statistics:'
      write(10,*)'program started ',time8()-istarttime,' s ago'
      write(10,*)'CPU time = ',int(etres), ' s (of which ',
     $     int(etres12(1)),' s user time)' 
      close(10)
      return
      end
********************************
      subroutine printvector(ip,p1,p2,qph)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3)
      common/ifirstprv/ifirst,icount
      data ifirst /1/
      data icount /0/
      icount = icount + 1
      if (ifirst.eq.1) then
         open(33,file='momenta',status='unknown')
         ifirst = 0
      else
         open(33,file='momenta',status='unknown',access='append')
      endif
      write(33,*)'point, count',ip,icount
      write(33,*)'set arrow 1 from 0,0,0 to',p1(3),',',p1(1),',',p1(2),
     .' lt 1'
      write(33,*)'set arrow 2 from 0,0,0 to',p2(3),',',p2(1),',',p2(2),
     .' lt 2'
      n = 0
      do k = 1,40
         if (qph(k,0).gt.0.d0) n = n+1
      enddo
      if (n.gt.0) then
         do k = 3,n+2
            q(0) = qph(k-2,0)
            q(1) = qph(k-2,1)
            q(2) = qph(k-2,2)
            q(3) = qph(k-2,3)
      write(33,*)'set arrow',k
     .,' from 0,0,0 to',q(3),',',q(1),',',q(2),''//
     .' lt 3'
        enddo
      endif
      close(33)
      return
      end
********************************
      subroutine mirror(p,q)
      double precision p(0:3),q(0:3)
      q(0) =  p(0)
      q(1) = -p(1)
      q(2) = -p(2)
      q(3) = -p(3)
      return
      end
********
      subroutine countphotons(q,eg,eta,na,nb)
      implicit double precision (a-h,o-z)
      dimension q(40,0:3),tmp(0:3)
      na = 0
      nb = 0
      do k = 1,40
         if (q(k,0).gt.eg) then
            na = na + 1
            do i = 0,3
               tmp(i) = q(k,i)
            enddo
            if (abs(rapidity(tmp)).lt.eta) nb = nb + 1
         endif
      enddo
      return
      end
**********
      subroutine new_boost(p,q,qq,idir)
      integer idir
      double precision p(0:3),q(0:3),qq(0:3),vboost(0:3)
      double precision b,g
! idir =  1: q is boosted where p is at rest and put in qq     
! idir = -1: q, in the frame where p is at rest, is boosted where p is not
!            at rest. The boosted q is put in qq.
      if (idir.eq.1) then
         vboost(0) =  0.d0
         vboost(1) =  p(1)/p(0)
         vboost(2) =  p(2)/p(0)
         vboost(3) =  p(3)/p(0)
      else
         vboost(0) =  0.d0
         vboost(1) = -p(1)/p(0)
         vboost(2) = -p(2)/p(0)
         vboost(3) = -p(3)/p(0)
      endif
      b=sqrt(vboost(1)**2+vboost(2)**2+vboost(3)**2)
      g=1.d0/sqrt(1.d0-b**2)
      call boost(g,vboost,q,qq)
      return
      end
*********
      double precision function lambda(x,y,z)
      double precision x,y,z
      lambda = x*x+y*y+z*z - 2.d0*x*y - 2.d0*x*z - 2.d0*y*z
      return
      end
*********
      double precision function theta(x)
      double precision x
      theta = 0.d0
      if (x.gt.0.d0) theta = 1.d0
      return
      end
**********
      function rapidity(p)
      implicit none
      double precision rapidity,p(0:3),p0mp3,p0pp3,ap,am
      p0mp3 = p(0)-p(3)
      p0pp3 = p(0)+p(3)
      ap = abs(p0pp3)
      am = abs(p0mp3)
      if (ap.le.0.d0) ap = 1.d-14
      if (am.le.0.d0) am = 1.d-14
      rapidity = 0.5d0*log(ap/am)
      return
      end
**********
      function pseudorapidity(p)
      implicit none
      double precision pseudorapidity,p(0:3),theta,tridot
      external tridot
      theta = acos(p(3)/sqrt(tridot(p,p)))
      pseudorapidity = -log(tan(theta/2.d0))
      return
      end
*********
      function nfactorial(n)
      integer n,nl,k,nfactorial      
      nl = n
      nfactorial = 1
      if (nl.gt.12) then
         print*,'factorial too big!!'
         return
      endif
      do while(nl.gt.0)
         nfactorial = nl * nfactorial
         nl = nl - 1
      enddo
      return
      end
**********
      function realpart(c)
      double precision realpart
      double complex c
      realpart = c
      return
      end
**********
      function e_(q1,q2,q3,q4)
      implicit double precision (a-h,o-z)
      dimension q1(0:3),q2(0:3),q3(0:3),q4(0:3)
      complex*16 esign,e_
      common/esig/esign
      esign = -(0.d0,1.d0)
c      e_ = 0.d0
c      return
!!! va moltiplicata per -1 rispetto a Form!!!
      e_ =q1(1)*q2(2)*q3(3)*q4(0)-q1(1)*q2(2)*q3(0)*q4(3)+
     .    q1(1)*q2(3)*q3(0)*q4(2)-q1(1)*q2(3)*q3(2)*q4(0)+
     .    q1(1)*q2(0)*q3(2)*q4(3)-q1(1)*q2(0)*q3(3)*q4(2)+
     .    q1(2)*q2(1)*q3(0)*q4(3)-q1(2)*q2(1)*q3(3)*q4(0)+
     .    q1(2)*q2(3)*q3(1)*q4(0)-q1(2)*q2(3)*q3(0)*q4(1)+
     .    q1(2)*q2(0)*q3(3)*q4(1)-q1(2)*q2(0)*q3(1)*q4(3)+
     .    q1(3)*q2(1)*q3(2)*q4(0)-q1(3)*q2(1)*q3(0)*q4(2)+
     .    q1(3)*q2(2)*q3(0)*q4(1)-q1(3)*q2(2)*q3(1)*q4(0)+
     .    q1(3)*q2(0)*q3(1)*q4(2)-q1(3)*q2(0)*q3(2)*q4(1)+
     .    q1(0)*q2(1)*q3(3)*q4(2)-q1(0)*q2(1)*q3(2)*q4(3)+
     .    q1(0)*q2(2)*q3(1)*q4(3)-q1(0)*q2(2)*q3(3)*q4(1)+
     .    q1(0)*q2(3)*q3(2)*q4(1)-q1(0)*q2(3)*q3(1)*q4(2)
       e_ = esign * e_ !!!! Carlo
      return
      end
*****************************************************
      subroutine messages(i)
      implicit double precision (a-h,o-z)
      parameter (nmesg = 8)  ! to be changed when adding a message!
      character w*(*),dp*(*)
      common/message_first/mesg_first(nmesg)
      data mesg_first /nmesg * 1/
      parameter (w='WARNING (')
      parameter (dp='): ')
      if (i.eq.1.and.mesg_first(1).eq.1) then
         print*,w,i,dp//'kinematics not suited for IS radiation!'
         mesg_first(1) = 0
      endif
      if (i.eq.2.and.mesg_first(2).eq.1) then
         print*,w,i,dp//'DO NOT USE IT FOR Z PRODUCTION!'
         mesg_first(2) = 0
      endif
      if (i.eq.3.and.mesg_first(3).eq.1) then
         print*,w,i,dp//'wrong choice for ALPHA!'
         mesg_first(3) = 0
      endif
      if (i.eq.4.and.mesg_first(4).eq.1) then
         print*,w,i,dp//'"hit or miss" warning:'
         print*,'=================='
c         print*,'i_over > 0: unweighted events are "not good"! '
c         print*,'Consider only weighted events/distributions'
         print*,'i_over > 0: some unweightening bias will be present!'
         print*,'=================='
         mesg_first(4) = 0
      endif
      if (i.eq.5.and.mesg_first(5).eq.1) then
         print*,' '
         print*,w,i,dp//'Gamma_W is fixed!!'
         print*,' '
         mesg_first(5) = 0
      endif
      if (i.eq.6.and.mesg_first(6).eq.1) then
         print*,w,i,dp//'in kinematics.f: IS RADIATION?!?!?'
         mesg_first(6) = 0
      endif
      if (i.eq.7.and.mesg_first(7).eq.1) then
         print*,' '
         print*,w,i,dp//'FIXED WIDTH!'
         mesg_first(7) = 0
         print*,' '
      endif
      if (i.eq.8.and.mesg_first(8).eq.1) then
         print*,' '
         print*,' '
         print*,w,i,dp//'RUNNING WIDTH!'
         mesg_first(8) = 0
         print*,' '
         print*,' '
      endif
!      if (i.eq. .and.mesg_first( ).eq.1) then
!         print*,w,i,dp//''
!         mesg_first( ) = 0
!      endif
      return
      end
******************************************************
      function cstar(p1,p2)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),psum(0:3)
      integer hadr1,hadr2
      common/is/hadr1,hadr2,i_pdf,iquark1,iquark2
! I need this common because for pp or ppbar cstar is differently defined
      do k=0,3
         psum(k) = p1(k) + p2(k)
      enddo
      rq2 = sqrt(2.d0)
! Collins - Soper momenta for particle 1 and 2 
      cs1p = (p1(0) + p1(3))/rq2
      cs2p = (p2(0) + p2(3))/rq2
      cs1m = (p1(0) - p1(3))/rq2
      cs2m = (p2(0) - p2(3))/rq2
      qmass = sqrt(dot(psum,psum))
      pt2 = psum(1)**2 + psum(2)**2
      cstar = 2.d0/qmass/sqrt(qmass**2 + pt2)*(cs1p*cs2m - cs1m*cs2p)
      if (hadr1.eq.hadr2) then
         sig = 1.d0
         if (psum(3).ne.0.d0) sig = abs(psum(3))/psum(3)
         cstar = cstar * sig
      endif
      return
      end
*-----------------------------------------------
       subroutine boost(g,v,p,q)
       implicit double precision (a-h,o-z)       
       double precision p(0:3),q(0:3),pp(0:3),v(0:3)
       do i=0,3
          pp(i)=p(i)
          q(i)=p(i)
       enddo
       if (g.eq.1.d0) return
       ppdv=tridot(pp,v)
       v2=v(1)**2+v(2)**2+v(3)**2
       q(0)=g*(pp(0)-ppdv)
       do i=1,3
          q(i)=pp(i)+(g-1.d0)*ppdv/v2*v(i)-g*v(i)*pp(0)
       enddo
       end
*------------------------------------------------
       function dot(p,q)
       implicit double precision (a-h,o-z)
       double precision p(0:3),q(0:3)       
       dot=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)       
       return
       end
*-----------------------------------------------
       function tridot(p,q)
       implicit double precision (a-h,o-z)
       double precision p(0:3),q(0:3)       
       tridot=p(1)*q(1)+p(2)*q(2)+p(3)*q(3)       
       return
       end
*-----------------------------------------------
      subroutine exchange_mom(p1,p2)
      implicit double precision (A-H,O-Z)
      dimension p1(0:3),p2(0:3),p4(0:3)
      do i=0,3
         p4(i)=p1(i)
      enddo
      do i=0,3
         p1(i)=p2(i)
      enddo
      do i=0,3
         p2(i)=p4(i)
      enddo      
      return
      end
******************************************************************
      subroutine rot(idir,vect,pin,pout)
      implicit double precision (a-h,o-z)       
      double precision pin(0:3),pout(0:3),pp(0:3),r(3,3),
     >                 vers(3),vect(0:3)
* This subroutine rotates the 4-vector pin in the frame where the z-axis is
* directed along the 4-vector vect(0,1,2,3). The rotated vector is stored
* in pout
* idir =  1 ---> direct rotation matrix
* idir = -1 ---> inverse rotation matrix
      pp(0) = pin(0)
      pp(1) = pin(1)
      pp(2) = pin(2)
      pp(3) = pin(3)
      vmo = 1.d0/sqrt(vect(1)**2+vect(2)**2+vect(3)**2)
      vers(1) = vect(1)*vmo
      vers(2) = vect(2)*vmo
      vers(3) = vect(3)*vmo
      vt = sqrt(vers(1)**2+vers(2)**2)
!   BUG - pointed out by CLEO people
!      v1ovt = vers(1)/vt
!      if (vt.eq.0.d0) v1ovt = 0.d0
!      v2ovt = vers(2)/vt
!      if (vt.eq.0.d0) v2ovt = 1.d0
      v1ovt = 0.d0
      v2ovt = 1.d0
      if (vt.gt.0.d0) then
         v1ovt = vers(1)/vt
         v2ovt = vers(2)/vt
      endif
      if (idir.eq.(-1)) then !! INVERSE rotation matrix
         r(1,1) =  vers(3)*v1ovt
         r(1,2) = -v2ovt
         r(1,3) =  vers(1)
         r(2,1) =  vers(3)*v2ovt
         r(2,2) =  v1ovt
         r(2,3) =  vers(2)
         r(3,1) = -vt
         r(3,2) =  0.d0
         r(3,3) =  vers(3)      
      else ! if (idir.eq.1) !! DIRECT rotation matrix
         r(1,1) =  vers(3)*v1ovt
         r(2,1) = -v2ovt
         r(3,1) =  vers(1)
         r(1,2) =  vers(3)*v2ovt
         r(2,2) =  v1ovt
         r(3,2) =  vers(2)
         r(1,3) = -vt
         r(2,3) =  0.d0
         r(3,3) =  vers(3)      
      endif
      pout(0) = pp(0)
      pout(1) = r(1,1)*pp(1) + r(1,2)*pp(2) + r(1,3)*pp(3)
      pout(2) = r(2,1)*pp(1) + r(2,2)*pp(2) + r(2,3)*pp(3)
      pout(3) = r(3,1)*pp(1) + r(3,2)*pp(2) + r(3,3)*pp(3)
      return
      end
************************************
      subroutine pick_2g(q,q1,q2)
      implicit double precision (a-h,o-z)
      dimension q(40,0:3),q1(0:3),q2(0:3)
*  LEADING ENERGETIC PHOTON IS EXTRACTED
        do i = 0,3
           q1(i) = 0.d0
           q2(i) = 0.d0
        enddo
        sum = q(1,0)+q(11,0)+q(21,0)+q(31,0)
        if (sum.lt.1.d-11) return
        enphot = q(1,0)
        enhard = enphot
        j=1
        do i = 1,40
           enphot = q(i,0)
           if (enhard.ge.enphot) then
              enhard = enhard
           else
              j = i
              enhard = enphot
           endif
        enddo 
        enhardl = enhard
        jl = j
        if (j.gt.0) then
           do i = 0,3
              q1(i) = q(j,i)
           enddo
        endif
*  next-to-leading energetic photon is extracted
      enhard = 0.d0
      j = 0
      do i = 1,40
         if (i.ne.jl) then
            enphot = q(i,0)
            if (enhard.ge.enphot) then
               enhard = enhard
            else
               j = i
               enhard = enphot
            endif
         endif
      enddo 
      enhard = enhard
      if (j.gt.0) then
         do i = 0,3
            q2(i) = q(j,i)
         enddo
      endif
      return
      end
**************************************
      function getphi(ppp)
      include 'shared.inc'
      dimension ppp(0:3)
      pm = sqrt(ppp(1)**2+ppp(2)**2+ppp(3)**2)
      c  = ppp(3)/pm
      s  = sqrt(1.d0-c**2)
      if (s.eq.0.d0) then
         getphi = 0.d0
         return
      else   
         arg = ppp(1)/pm/s
*  avoiding numerical problems......
         if (abs(arg).ge.1.d0) then
            iarg = arg
            arg  = iarg
         endif   
         if (ppp(2).ge.0.d0) getphi = acos(arg)
         if (ppp(2).lt.0.d0) getphi = 2.d0*pi-acos(arg)
      endif             
      return
      end
***
      function deltaalpha(q2)
      include 'shared.inc'
      real*4 qin4,st24,der4,errder4,deg4,errdeg4
      dimension amasses(9)
* e,mu,tau,top masses and light quarks !!:
      amasses(1) = me
      amasses(2) = mmu
      amasses(3) = mtau
      amasses(4) = mt
      amasses(5) = mu
      amasses(6) = mc
      amasses(7) = md
      amasses(8) = ms
      amasses(9) = mb
      somma=0.d0
      do i=1,9
         somma=somma+summa(amasses(i),q2,i)
      enddo
      deltaalpha=alpha/pi*somma
      return
      end
*-----------------------------------------------------
* vacuum polarization: alpha(0)->alpha(0)*vpol(q2)
      function vpol(q2) 
      include 'shared.inc'
      real*4 qin4,st24,der4,errder4,deg4,errdeg4
      dimension amasses(4)
      common/vpolfirst/ifirst
      data ifirst/1/
* e,mu,tau,top masses: 
      amasses(1) = me
      amasses(2) = mmu
      amasses(3) = mtau
      amasses(4) = mt
*
      somma=0.d0
      do i=1,4
         somma=somma+summa(amasses(i),q2,i)
      enddo
      amodq2=abs(q2)
      st2=s2tw
      qin=sqrt(amodq2)
      qin4=qin
      st24=st2
      if (qin4.ge.1000.d0) then
         qin4 = 999.9999d0
      endif
      call hadr5n(qin4,st24,der4,errder4,deg4,errdeg4)
      der=der4
      deg=deg4
      errder=errder4
      errdeg=errdeg4
      dalphaqui=alpha/pi*somma+der
      vpol=1.d0/(1.d0-dalphaqui)
      ifirst = 0
      return
      end
*-----------------------------------------------------
* vacuum polarization: alpha(0)->alpha(0)*vpol(q2)
      function vpolan(q2) 
      include 'shared.inc'
      real*4 qin4,st24,der4,errder4,deg4,errdeg4
      dimension amasses(4)
      common/vpolnewfirst/ifirst
      data ifirst/1/
*
      if (ifirst.eq.1) then
         do k=1,10
c            print*,'WARNING IN VPOL: dalpha5 not from FJ!'
         enddo
      endif
      dalphaqui = deltaalpha(q2) !! dalpha is prop. to alpha!!
      vpolan   = 1.d0/(1.d0-dalphaqui)
      ifirst    = 0
      return
      end
*----------------------------------------------------------
* leptonic and top contribution to vacuum polarization
      function summa(am,q2,i)
      implicit double precision (a-h,o-z)
      dimension nc(9),qf2(9)
* nc and qf are color factor (1 for leptons, 
* 3 for quarks) and charge
      do j=1,3
         nc(j) = 1.d0
         qf2(j)= 1.d0
      enddo
      do j = 4,6
         nc(j) = 3.d0
         qf2(j)= 4.d0/9.d0
      enddo
      do j = 7,9
         nc(j) = 3.d0
         qf2(j)= 1.d0/9.d0
      enddo
      am2=am**2
      if (q2.ge.0.d0.and.q2.lt.(4.d0*am2)) then
         sq=sqrt(4.d0*am2/q2-1.d0)
         summa=nc(i)*qf2(i)*(-5.d0/9.d0-(4.d0/3.d0)*(am2/q2)+
     >        (4.d0/3.d0*(am2/q2)**2+1.d0/3.d0*am2/q2-1.d0/6.d0)*
     >        4.d0/sq*atan(1.d0/sq))
      else
         sq=sqrt(1.d0-4.d0*am2/q2)
         arglog=abs((1.d0-sq)/(1.d0+sq))
         summa=nc(i)*qf2(i)*(-5.d0/9.d0-(4.d0/3.d0)*(am2/q2)+
     >        (4.d0/3.d0*(am2/q2)**2+1.d0/3.d0*am2/q2-1.d0/6.d0)*
     >        2.d0/sq*log(arglog))
      endif
      return
      end
****************************************************
      REAL*8 FUNCTION DDILOGhere(X)
      REAL*8 X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      REAL*8 C(0:18),H,ALFA,B0,B1,B2
      DATA ZERO /0.0D0/, ONE /1.0D0/
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/
      DATA C( 0) / 0.42996 69356 08137 0D0/
      DATA C( 1) / 0.40975 98753 30771 1D0/
      DATA C( 2) /-0.01858 84366 50146 0D0/
      DATA C( 3) / 0.00145 75108 40622 7D0/
      DATA C( 4) /-0.00014 30418 44423 4D0/
      DATA C( 5) / 0.00001 58841 55418 8D0/
      DATA C( 6) /-0.00000 19078 49593 9D0/
      DATA C( 7) / 0.00000 02419 51808 5D0/
      DATA C( 8) /-0.00000 00319 33412 7D0/
      DATA C( 9) / 0.00000 00043 45450 6D0/
      DATA C(10) /-0.00000 00006 05784 8D0/
      DATA C(11) / 0.00000 00000 86121 0D0/
      DATA C(12) /-0.00000 00000 12443 3D0/
      DATA C(13) / 0.00000 00000 01822 6D0/
      DATA C(14) /-0.00000 00000 00270 1D0/
      DATA C(15) / 0.00000 00000 00040 4D0/
      DATA C(16) /-0.00000 00000 00006 1D0/
      DATA C(17) / 0.00000 00000 00000 9D0/
      DATA C(18) /-0.00000 00000 00000 1D0/
      IF(X .EQ. ONE) THEN
         DDILOGhere=PI6
         RETURN
      ELSE IF(X .EQ. MONE) THEN
         DDILOGhere=MALF*PI6
         RETURN
      END IF
      T=-X
      IF(T .LE. MTWO) THEN
         Y=MONE/(ONE+T)
         S=ONE
         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
      ELSE IF(T .LT. MONE) THEN
         Y=MONE-T
         S=MONE
         A=LOG(-T)
         A=-PI6+A*(A+LOG(ONE+ONE/T))
      ELSE IF(T .LE. MALF) THEN
         Y=(MONE-T)/T
         S=ONE
         A=LOG(-T)
         A=-PI6+A*(MALF*A+LOG(ONE+T))
      ELSE IF(T .LT. ZERO) THEN
         Y=-T/(ONE+T)
         S=MONE
         A=HALF*LOG(ONE+T)**2
      ELSE IF(T .LE. ONE) THEN
         Y=T
         S=ONE
         A=ZERO
      ELSE
         Y=ONE/T
         S=MONE
         A=PI6+HALF*LOG(T)**2
      END IF
      H=Y+Y-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 1 I = 18,0,-1
         B0=C(I)+ALFA*B1-B2
         B2=B1
    1 B1=B0
      DDILOGhere=-(S*(B0-H*B2)+A)
      RETURN
      END
