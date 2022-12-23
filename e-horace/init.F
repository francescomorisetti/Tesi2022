      subroutine userinterface
      include 'shared.inc'
      character*10 tmpchar
      character*6 ord
      common/qedORDER/ord
      common/nphot_mode/nphotmode
      common/reading/ihavetoread,iparallel
      common/toreadsubtraction/ireadsub
      character*20 readinput

      print*,' '
      print*,' *****************************************'
      print*,' *                                       *'
      print*,' *         Welcome to Horace ',horaceversion,'       *'
      print*,' *                                       *'
      print*,' *****************************************'
      print*,'http://www.pv.infn.it/hepcomplex/horace.html'
      print*,' '

      
****  USER INPUTS
c      print*,'Which pdf replica? (int 0-100)'
c      read*,ipdfreplica
      ipdfreplica = 0
      call resetname(cpdfreplica)
      call itoa(ipdfreplica,cpdfreplica)
      ilrep = 3
      if (ipdfreplica.lt.100) ilrep = 2
      if (ipdfreplica.lt.10)  ilrep = 1      
      
      print*,'Start new (and best) or old HORACE mode? (new/old)'
      read*,tmpchar
      iewk = 1
      if (tmpchar.eq.'old') iewk = 0
      print*,'Select the calculation of rad. corr. (born/alpha/exp)'
      read(*,'(A)')ord
      print*,'Enter initial state hadrons'
      print*,'(2212 = electron, -2212 = positron)'
      read*,hadr1
      read*,hadr2
      print*,'Select the PDF set according to HORACE scheme'
      read*,iwhich
      print*,'Enter the center of mass energy (GeV)'
      read*,ebeam
      ebeam = ebeam/2.d0
      print*,'Select Drell-Yan boson production (W+/W-/Z/W+andW-)'
      read*,tmpchar
      isok = 0
      if (tmpchar.eq.'W+andW-') then
         ichfs = 2
         isok  = 1
      elseif (tmpchar.eq.'W+') then
         ichfs = 1
         isok  = 1
      elseif (tmpchar.eq.'W-') then
         ichfs = -1
         isok  =  1
      elseif (tmpchar.eq.'Z') then
         ichfs = 0
         isok  = 1
      endif
      if (isok.eq.0) then 
         print*,'Select correctly the boson production!'
         stop
      endif
      chfs  = min(ichfs,1)
      boson = 'W'
      if (abs(chfs).lt.1.d0) boson = 'Z'

c      print*,'Enter the EW input scheme: 0=alpha(0), 1=Gmu, 2=IBA'
      print*,'Enter the EW input scheme (0/1)'
      print*,' 0 --> pure alpha(0) scheme'
      print*,' 1 --> Gmu scheme for W and Z production '
      print*,' 2 --> best (only for Z production)'
      read*,ischemealpha
      if (ischemealpha.eq.2.and.boson.eq.'W') ischemealpha = 1

      print*,'Select PDF scale choice (1/2/3)'
      print*,' 1 --> fixed to M_boson'
      print*,' 2 --> running, sqrt(M_leptons^2)'
      print*,' 3 --> running, sqrt(M_leptons^2 + pt_lepton_pair^2)'
      read*,iscalechoice
      
      print*,'Select decay mode (electron/muon/tau)'
      read*,tmpchar
      lepton = 1
      if (tmpchar.eq.'muon') lepton = 2
      if (tmpchar.eq.'tau') lepton = 3
      print*,'Enter the directory where to store the ouput files'
      read(*,'(A)')path
      print*,'Enter the number of event to generate'
      read*,nmax
      print*,'Events are weighted or unweighted? (weighted/unweighted)'
      read(*,'(A)')tmpchar
      idumpweighted = 1
      if (tmpchar.eq.'unweighted') idumpweighted = 0
      print*,'Enter the upper limit for unweightening'
      read*,tune_hm
      print*,'Events have to be stored? (yes/no)'
      read(*,'(A)')tmpchar
      storing = tmpchar(1:1)
      print*,'Calculate or read subtraction terms? (calculate/read) '
      read(*,'(A)')tmpchar
      ireadsub = 1
      if (tmpchar.eq.'calculate') ireadsub = 0
      print*,'External file to read or write subtraction terms'
      read*,subtractionfile
      print*,'Include photon-induced processes? (yes/no)'
      print*,'(effective only if using PDFs with photons)'
      read*,tmpchar
      iphinduced = 0
      if (tmpchar.eq.'yes') iphinduced = 1

      if (iewk.eq.0) iphinduced = 0

      print*,'Enter seed for the random number generator'
      read*,iseed
c      print*,'Input for cuts:'
c      print*,'Enter maximum lepton(s) rapidity'
c      read*,etamax
c      print*,'Enter minimum lepton pt (GeV)'
c      read*,ptmin
c      print*,'Enter minimum lepton (missing) pt (GeV)'
c      read*,ptminmiss
c      tmmin = 0.d0
c      print*,'Apply a cut on transverse mass (inv. mass for Z)? (0/1)'
c      read*,icutmt
c      print*,'Enter the minimum mt cut (active if previous = 1) (GeV)'
c      read*,tmmin

      
* Parameters
      pi      = 4.d0*atan(1.d0)
      convfac = 0.389379660D9
      
* Masses (GeV)
      me   =  0.510998928d-3 ! electron mass
      mmu  =  0.1056583715d0   ! muon mass
      mtau =  1.77682d0     ! tau mass

      mh   = 125.d0          ! Higgs mass
      mu   =   0.06983d0     ! up quark mass
      md   =   0.06984d0     ! down quark mass
      ms   =   0.150d0       ! strange quark mass
      mc   =   1.2d0         ! charm quark mass
      mb   =   4.6d0         ! bottom quark mass
      mt   = 173.5d0         ! top quark mass

* W and Z masses and total widths (GeV)

** REPORT TUNED COMPARISON (AND DEFAULT) SETUP
      mw = 80.385d0       ! W mass
      mz = 91.1876d0      ! Z mass
      gw =  2.085d0
      gz =  2.4952d0

#ifdef REPORT
** REPORT BENCHMARK SETUP
      if (ischemealpha.eq.2) ischemealpha = 1
      if (ischemealpha.eq.1) then
         mw = 80.358d0
         mz = 91.1535d0
         gw =  2.084d0
         gz =  2.4943d0
      endif
#endif
****  

*  Constants
      alpha      = 1.d0/137.035999074d0   ! alpha_em
      gf         = 1.1663787d-5           ! GeV^(-2) ! Fermi constant
      s2tw       = 1.d0 - mw*mw/mz/mz     ! sin^2 of the weak mixing angle
      s2twlepexp = 0.23152d0              ! sin^2 t_eff_leptonic
      alpha_gf   = sqrt(2.d0)*gf*mw*mw*s2tw/pi

*************************************************
* CKM matrix elements (only the absolute value of the matrix elements)
*   V_ud V_us V_ub
*   V_cd V_cs V_cb
*   V_td V_ts V_tb
      ckm(1,1) = 0.975d0 ! V_ud
      ckm(2,1) = 0.222d0 ! V_cd
      ckm(3,1) = 0.d0    ! V_td
      ckm(1,2) = 0.222d0 ! V_us
      ckm(2,2) = 0.975d0 ! V_cs
      ckm(3,2) = 0.d0    ! V_ts
      ckm(1,3) = 0.d0    ! V_ub
      ckm(2,3) = 0.d0    ! V_cb
      ckm(3,3) = 1.d0    ! V_tb
*******************************************************
      return
      end
*****************************************
      subroutine init
* This soubroutine initialize various variables and puts them in common blocks
      include 'shared.inc'
      character*10 model,strwunw
      character*6 ord
      character*7 strprod,strfs
      common/reading/ihavetoread,iparallel
      common/deltaerre/deltar
      character*30 whichdeltar
      common/which_dr/whichdeltar
      common/nphot_mode/nphotmode
      common/input_matching/model
      common/qedORDER/ord
      common/regulators/ESOFTMAX
      common/reject/iswitchreject
      common/switchdumperrors/idumperrors
      double complex mwcomplex
      common/resonant/mwcomplex

#include "looptools.h"
      
*[[[  HORACE VERSION
      horaceversion = '3.2'
*]]]
      do k = 1,100
         outfile(k:k)  = ' '
         noefile(k:k)  = ' '
         path(k:k)     = ' '
         cutsfile(k:k) = ' '
      enddo
      do k = 1,150
         storfile(k:k) = ' '
         conffile(k:k) = ' '
         command(k:k)  = ' '
      enddo
************************************************************
      if (ihavetoread.eq.1) call userinterface
************************************************************

**** once these were in input...
      ismear = 0
      irec   = 0
c      icutmt = 0
c      tmmin  = 0.d0
****
      if (lepton.eq.1) mfs = me 
      if (lepton.eq.2) mfs = mmu
      if (lepton.eq.3) mfs = mtau
*** Inner input parameters, which the user may not need to modify
      iwriteoutput =  10000
      idumperrors  =  0

      nphotmode = -1
      
      i_pdf   = 1 ! calculate hadron level cross section
      iquark1 = 0
      iquark2 = 0
      iplsub  = 1

c$$$      print*,'PARTON LEVEL WITHOUT SUBTRACTION!!!'
c$$$      i_pdf   = 0 ! calculate hadron level cross section
c$$$      iquark1 = 2
c$$$      iquark2 = 0
c$$$      iplsub = 0
c$$$      nphotmode = 1
***
      iswitchreject = 0
      
      iswitchreject = 1
      
      qedorder      = 'exp'
      if (ord.eq.'alpha') qedorder = 'alpha'
      ptqed    = 'yes'
      model    = 'matchedps'

      rescalepdfscale = 1.d0

****  
      charge_factor(1) = 0
      charge_factor(2) = 0
      charge_factor(3) = 1
      if (ord.eq.'born') charge_factor(3) = 0
      charge_factor(4) = 0
      if (boson.eq.'Z')  charge_factor(4) = 1
      if (ord.eq.'born') charge_factor(4) = 0
**** QED Parton Shower parameters...
      eps = 4.d-4
      if (boson.eq.'Z') eps=7.d-4
      if (qedorder.eq.'exp') eps = 1.d-7
** ALPHA for Z production
      ialphachoice = 2
      ialphachoice = 0
*--- the minimum energy of the photon is eps*ebeam (and not eps*ecms!)
      eps = 1.d-5
      if (boson.eq.'Z'.and.qedorder.eq.'alpha') eps=1.d-3

*** FIXED
      esoftmax = 0.004d0
      if (i_pdf.eq.0) esoftmax = eps*ebeam

      i_scheme = 4
      new_kine = 1
**************************************************
      iplus  = -2.d0*log(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
      ebeam1 = ebeam
      ebeam2 = ebeam
*------------------------
* For conventions, I set (PDG MC numbering scheme):
*   1 = d
*   2 = u
*   3 = s
*   4 = c
*   5 = b
* - 1 = d bar
* - 2 = u bar
* - 3 = s bar
* - 4 = c bar
* - 5 = b bar
* In mq(i) and chq(i), I load the masses and charges of the corresponding quark
      mq(1)  = md
      mq(2)  = mu
      mq(3)  = ms
      mq(4)  = mc
      mq(5)  = mb
      mq(0)  = 0.d0
      mq(-1) = md
      mq(-2) = mu
      mq(-3) = ms
      mq(-4) = mc
      mq(-5) = mb
cccccccccccccccccccccccccccccccccccccc
c**      mq(i) all set to mu or md!!!      
      mq(1)  = md
      mq(2)  = mu
      mq(3)  = md
      mq(4)  = mu
      mq(5)  = md
      mq(0)  = 0.d0
      mq(-1) = md
      mq(-2) = mu
      mq(-3) = md
      mq(-4) = mu
      mq(-5) = mu
cccccccccccccccccccccccccccccccccccccc
      chq(1)  =  -1.d0/3.d0
      chq(2)  =   2.d0/3.d0
      chq(3)  =  -1.d0/3.d0
      chq(4)  =   2.d0/3.d0
      chq(5)  =  -1.d0/3.d0
      chq(0)  =   0.d0
      chq(-1) =   1.d0/3.d0
      chq(-2) =  -2.d0/3.d0
      chq(-3) =   1.d0/3.d0
      chq(-4) =  -2.d0/3.d0
      chq(-5) =   1.d0/3.d0
      itre(1)  = -1.d0/2.d0
      itre(2)  =  1.d0/2.d0
      itre(3)  = -1.d0/2.d0
      itre(4)  =  1.d0/2.d0
      itre(5)  = -1.d0/2.d0
      itre(0)  =  0.d0
      itre(-1) =  1.d0/2.d0
      itre(-2) = -1.d0/2.d0
      itre(-3) =  1.d0/2.d0
      itre(-4) = -1.d0/2.d0
      itre(-5) =  1.d0/2.d0
      do j = -ifl, ifl
         do i = -ifl, ifl
            myckm(i,j) = 0.d0
         enddo
      enddo
!d s b | 1 3 5
!u c t | 2 4 6
      myckm(2,-1)  = ckm(1,1) !V_ud
      myckm(1,-2)  = ckm(1,1)
      myckm(-1,2)  = ckm(1,1)
      myckm(-2,1)  = ckm(1,1)
      myckm(1,2)   = ckm(1,1)
      myckm(2,1)   = ckm(1,1)
      myckm(-1,-2) = ckm(1,1)
      myckm(-2,-1) = ckm(1,1)

      myckm(2,-3)  = ckm(1,2) !V_us
      myckm(3,-2)  = ckm(1,2) 
      myckm(-2,3)  = ckm(1,2)
      myckm(-3,2)  = ckm(1,2)
      myckm(2,3)   = ckm(1,2)
      myckm(3,2)   = ckm(1,2)
      myckm(-2,-3) = ckm(1,2)
      myckm(-3,-2) = ckm(1,2)

      myckm(2,-5)  = ckm(1,3) !V_ub
      myckm(5,-2)  = ckm(1,3)
      myckm(-2,5)  = ckm(1,3)
      myckm(-5,2)  = ckm(1,3)
      myckm(2,5)   = ckm(1,3)
      myckm(5,2)   = ckm(1,3)
      myckm(-2,-5) = ckm(1,3)
      myckm(-5,-2) = ckm(1,3)

      myckm(4,-1)  = ckm(2,1) !V_cd
      myckm(1,-4)  = ckm(2,1)
      myckm(-4,1)  = ckm(2,1)
      myckm(-1,4)  = ckm(2,1)
      myckm(1,4)   = ckm(2,1)
      myckm(4,1)   = ckm(2,1)
      myckm(-1,-4) = ckm(2,1)
      myckm(-4,-1) = ckm(2,1)

      myckm(4,-3)  = ckm(2,2) !V_cs
      myckm(3,-4)  = ckm(2,2)
      myckm(-4,3)  = ckm(2,2)
      myckm(-3,4)  = ckm(2,2)
      myckm(3,4)   = ckm(2,2)
      myckm(4,3)   = ckm(2,2)
      myckm(-3,-4) = ckm(2,2)
      myckm(-4,-3) = ckm(2,2)

      myckm(4,-5)  = ckm(2,3) !V_cb
      myckm(5,-4)  = ckm(2,3)
      myckm(-4,5)  = ckm(2,3)
      myckm(-5,4)  = ckm(2,3)
      myckm(4,5)   = ckm(2,3)
      myckm(5,4)   = ckm(2,3)
      myckm(-4,-5) = ckm(2,3)
      myckm(-5,-4) = ckm(2,3)

** "CKM" for photon induced
      myckm(0,0)  = 0.d0
      myckm(0,1)  = sqrt(myckm(1,2)**2+myckm(1,4)**2)
      myckm(0,2)  = sqrt(myckm(2,1)**2+myckm(2,3)**2+myckm(2,5)**2)
      myckm(0,3)  = sqrt(myckm(3,2)**2+myckm(3,4)**2)
      myckm(0,4)  = sqrt(myckm(4,1)**2+myckm(4,3)**2+myckm(4,5)**2)
      myckm(0,5)  = sqrt(myckm(5,2)**2+myckm(5,4)**2)
      myckm(0,-1) = sqrt(myckm(1,2)**2+myckm(1,4)**2)
      myckm(0,-2) = sqrt(myckm(2,1)**2+myckm(2,3)**2+myckm(2,5)**2)
      myckm(0,-3) = sqrt(myckm(3,2)**2+myckm(3,4)**2)
      myckm(0,-4) = sqrt(myckm(4,1)**2+myckm(4,3)**2+myckm(4,5)**2)
      myckm(0,-5) = sqrt(myckm(5,2)**2+myckm(5,4)**2)
      myckm(1,0)  = myckm(0,1)
      myckm(2,0)  = myckm(0,2)
      myckm(3,0)  = myckm(0,3)
      myckm(4,0)  = myckm(0,4)
      myckm(5,0)  = myckm(0,5)
      myckm(-1,0) = myckm(1,0)
      myckm(-2,0) = myckm(2,0)
      myckm(-3,0) = myckm(3,0)
      myckm(-4,0) = myckm(4,0)
      myckm(-5,0) = myckm(5,0)
*************************************

*******************************************************
* The random number generator is initialized
      call initrng(iseed)
* PDFs are initialized...
      call initPDFhorace
*******************************************************
* Initialization of various quantities in common blocks
      do j = -ifl, ifl
         do i = -ifl, ifl
            pdfpdf(i,j) = 0.d0
         enddo
         pdf1(j) = 0.d0
         pdf2(j) = 0.d0
      enddo 
      do k = 1,100
         if (path(k:k).ne.' ') lpath=k
      enddo
      if (path(lpath:lpath).ne.'/') then
         lpath = lpath+1
         path(lpath:lpath)='/' 
      endif
      open(11,file=path,status='old',iostat=ierr)
      if (ierr.eq.0) then
         print*,' '
         print*,'WARNING:' 
         print*,'The directory ',path(1:lpath),' already exists.'
         print*,'Overwriting the files in it...'
         print*,' '
      endif
      close(11)
      command(1:9) = 'mkdir -p '
      command(10:) = path
      call system(command)
      outfile(1:lpath)  = path(1:lpath) 
      noefile(1:lpath)  = path(1:lpath) 
      storfile(1:lpath) = path(1:lpath) 
      cutsfile(1:lpath) = path(1:lpath) 
      outfile(lpath+1:)  = 'statistics.txt' 
      noefile(lpath+1:)  = 'numberofevents.txt'
      cutsfile(lpath+1:) = 'cuts_and_fun.f'
c      if (storing.ne.'y') then 
c         storfile(lpath+1:) = 'empty.evts'
c         conffile(lpath+1:) = 'empty.conf'
c      else
      ip  = lnblnk(path)
      io  = lnblnk(ord)
      ivr = lnblnk(horaceversion)
      if (ichfs.eq. 2) strprod = 'W+W-'
      if (ichfs.eq. 1) strprod = 'W+'
      if (ichfs.eq.-1) strprod = 'W-'
      if (ichfs.eq. 0) strprod = 'Z'
      if (iphinduced.eq.1) then
         if (ichfs.eq. 2) strprod = 'W+W-_pi'
         if (ichfs.eq. 1) strprod = 'W+_pi'
         if (ichfs.eq.-1) strprod = 'W-_pi'
         if (ichfs.eq. 0) strprod = 'Z_pi'
      endif
      isp = lnblnk(strprod)
      if (lepton.eq.1) strfs = 'e'
      if (lepton.eq.2) strfs = 'mu'
      if (lepton.eq.3) strfs = 'tau'
      isf = lnblnk(strfs)
      if (idumpweighted.eq.1) then 
         strwunw = 'weighted'
      else
         strwunw = 'unweighted'
      endif
      iwunw = lnblnk(strwunw)
      storfile = path(1:ip)//'horace_'//horaceversion(1:ivr)
     $     //'_'//strprod(1:isp)//'_'//strfs(1:isf)//'_'
     $     //ord(1:io)//'_'//strwunw(1:iwunw)
c      conffile =
c     .  storfile(1:lnblnk(storfile))//'-'//cpdfreplica(1:ilrep)//'.conf'
c      outfile  =
c     .  storfile(1:lnblnk(storfile))//'-'//cpdfreplica(1:ilrep)//'.stat'
c      storfile =
c     .  storfile(1:lnblnk(storfile))//'-'//cpdfreplica(1:ilrep)//'.evts'
      conffile =
     .     storfile(1:lnblnk(storfile))//'.conf'
      outfile  =
     .  storfile(1:lnblnk(storfile))//'.stat'
      storfile =
     .  storfile(1:lnblnk(storfile))//'.evts'

!     creating symbolinc link for horace.conf
      if (storing.eq.'y') then
         icnf = lnblnk(conffile)
         command = 'ln -s -f '//conffile(ip+1:icnf)//' '//
     .        path(1:ip)//'horace.conf'
         call system(command)
      endif

      open(11,file=noefile,status='unknown')
      write(11,*)nmax
      close(11)
************* FEYNARTS INTERFACE & Dr
c      deltar = 0.028455588d0 ! from Bardin and Arbuzov
c      deltar = deltar * alpha/alpha_gf
*** needed to calculate deltar_r
      call fillfeynartscommon

      mwcomplex = dcmplx(mw,-gw/2.d0)
      if (boson.eq.'Z') mwcomplex = dcmplx(mz,-gz/2.d0)

      
c      call ffini
      call ltini

      call setversionkey(1*KeyD0)
c      call setversionkey(1*KeyAll)
      
      call setlambda(1.d-10)
      call setdelta(0.d0)
      call setmudim(1.d0)
      call calc_renconst   ! for W
      call uucalc_renconst ! for Z
      call ddcalc_renconst ! for Z
*********
      whichdeltar = ' (one-loop EW delta_r)'
      deltar = getdeltar()
***

      if (ischemealpha.lt.1) then
         deltar   = 0.d0
         alpha_gf = alpha
      endif

      if (storing.eq.'y') call initstorage
      
      return
      end
