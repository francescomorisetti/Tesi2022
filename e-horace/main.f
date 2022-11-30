      program HORACE
      include 'shared.inc'
      integer isvec(25)
      common/rlxisvec/isvec
      common/switchdumperrors/idumperrors
      real*4 rnd(1)
      integer*8 ncalls,ncut,nweighted,nevunwneg
      dimension x(40),xvect(10),vboost(0:3)
      common/reading/ihavetoread,iparallel
      common/reject/iswitchreject
      common/fordebugging/idebugging
      common/nevents/nweighted,nevunwneg
      common/n_phot/nphot
      common/radpattern/nph(4)
      common/fordumping/varw,varbefore,iev,idestroyed
      common/mcsamplingweight/wsmc
      common/tetto/tetto
      integer*8 istarttime
      common/starttime/istarttime
      data ((qph(i,j),i=1,40),j=0,3) /160 * 0.d0/
      istarttime = time8()
      idebugging  = 0
      ihavetoread = 1
      iparallel   = 0
*****************************
      call init
*****************************
      i_over    = 0
      iwrite    = 0
      ncalls    = 0
      nweighted = 0
      nevunwneg = 0
      wsum      = 0.d0 ! weighted sum for integrator without HIT or MISS
      wsum2     = 0.d0
      wsumneg   = 0.d0
      wsumneg2  = 0.d0
      wsectneg  = 0.d0
      wvarneg2  = 0.d0
      wsumov    = 0.d0
      wsumov2   = 0.d0
      wsectov   = 0.d0
      wvarov2   = 0.d0
      nneg      = 0
      xmax      = xpdfmax
      xmin      = xpdfmin
      
      if (i_pdf.eq.0) then
         xmax = 1.d0
         xmin = 0.d0
      endif
* The upper limit (tetto) for hit or miss is (very) roughtly guessed....
      amw2 =  mw**2
      c0   =  0.d0
      i1   =  1
      i2   =  2
      E1cm = mw / 2.d0
      E2cm = mw / 2.d0
      pz1  = E1cm
      pz2  = -E2cm
      ss   = 80.d0 !2.d0*ebeam1
      xx   = (mw-gw)/ss
      if (boson.eq.'Z') amw2 = mz**2
      if (boson.eq.'Z') xx =(mz-gz)/ss 
      if (xx.gt.1) xx = 0.8d0
      aln = 1.d0
      if (xmin.gt.0.d0) aln=log(xmax/xmin)

cc dummy call to cuts in order to initialise some quantities in common!
      ptmp(0) = 1.d0
      ptmp(1) =-0.001d0
      ptmp(2) = 0.001d0
      ptmp(3) = 0.99d0
      call cuts(ptmp,ptmp,qph,icutt)
ccc
      call get_pdf(xx,xx,ss)
      sup_pdf = pdf1(2)*xx*aln*1.5d0
      if (sup_pdf.lt.0.001.or.i_pdf.eq.0) sup_pdf = 1.d0
      tetto = shard(amw2,c0,2,-1)*4.d0*pi/4.d0*8.d0/3.d0
      fatt = 1.d0
      fbw = gw*mw/ss
      tetto = tetto*sup_pdf*sup_pdf*fatt*fbw*tune_hm
      if (boson.eq.'Z') tetto = tetto * 100.d0
      ratio_max = 0.d0
*=DEBUG=============================================
c$$$      open(16,file='ranluxsequence',status='unknown')
c$$$      read(16,*)isvec
c$$$      call rluxin(isvec)
c$$$      close(16)
*=DEBUG=============================================
      varbefore = 0.d0
      iev = 0
      ncut = 0
      istop = 0
* DO loop over events is started
      do while(istop.eq.0)
*=DEBUG=============================================
c$$$         if (iev.gt.2552115) then
c$$$            idebugging = 1
c$$$         endif
c$$$         iiiii = iev
c$$$         if (mod(iiiii,2000).eq.0) then
c$$$            print*,'@@@@@@@@@@@@@@@@@@@@@@@@@'
c$$$            print*,'======================================'
c$$$            print*,'DUMPING RANLUX SEQUENCE'
c$$$            print*,'======================================'
c$$$            print*,'@@@@@@@@@@@@@@@@@@@@@@@@@'
c$$$c            print*,'evento:',iev
c$$$            open(16,file='ranluxsequence',status='unknown')
c$$$            call rluxut(isvec)
c$$$            write(16,*)isvec
c$$$            close(16)
c$$$         endif
*=DEBUG=============================================
         i_reg = 1
         do while(i_reg.eq.1) ! do LOOP for HIT or MISS and exp cuts
            ncalls = ncalls+1

ccccccccccccc
            call rluxut(isvec)
ccccccccccccc
            
            if (idumpweighted.eq.1.and.nweighted.eq.nmax) then 
               istop = 1
               i_reg = 0
            endif
	    if (iewk.eq.0) then
	       call generate_event_OLDPS(sdif,i_ok,p3,p4,qph)
	    else               
               call generate_event_matched(sdif,i_ok,p3,p4,qph)
            endif
            if (abs(sdif).lt.1.d-30.or.i_ok.eq.0) ncut = ncut + 1
            if (i_ok.eq.0) then
               sdif = 0.d0
            endif

c            if (isnan(sdif)) then
c               print*,'main.f line 131'
c               print*,sdif,ncalls
c            endif
               
            if (i_ok.eq.1) then
               nweighted = nweighted + 1
               if (mod(nweighted,iwriteoutput).le.1d-5) iwrite = 1
               if (storing.eq.'y'.and.idumpweighted.eq.1) 
     .              call storehoraceevent(sdif,p3,p4,qph)
** The binning of the distributions is done here.....
               call distributions(sdif,ncalls,p3,p4,qph)
cc               call distributions2(sdif,ncalls,p3,p4,qph)
            endif
* The weighted sum is calculated
* Hit or miss for unweightening is done after
            wsum  = wsum + sdif
            wsum2 = wsum2+sdif**2
            wsect = wsum/ncalls
            wvar2before = wvar2
            wvar2 = (wsum2/ncalls-wsect**2)/ncalls
            if (i_ok.eq.1.and.sdif.lt.0.d0) then 
               nneg = nneg + 1
               wsumneg  = wsumneg + sdif
               wsectneg = wsumneg/ncalls 
               wsumneg2 = wsumneg2 + sdif**2
               wvarneg2 = (wsumneg2/ncalls - wsectneg**2)/ncalls
            endif
            varw = sqrt(abs(wvar2))

            if (varw.gt.1.1d0*varbefore) then
               idestroyed = 0
c               print*,'===== VARIANCE BUMP ======'
c               print*,varw,varbefore,iev,nphot
               if (ncalls.gt.10000000.and.iswitchreject.eq.1
     .              .and.varw.gt.2d0*varbefore) then
                  print*,'>>> REJECTING THE EVENT!! (PID=',getpid(),')'
                  idestroyed = 1
                  ncalls = ncalls - 1
                  iev    = iev - 1
                  wsum   = wsum - sdif
                  wsum2  = wsum2 - sdif**2
                  wsect  = wsum/ncalls
                  wvar2  = (wsum2/ncalls-wsect**2)/ncalls
                  if (sdif.lt.0.d0) then
                     nneg = nneg -1
                     wsumneg = wsumneg - sdif
                     wsumneg2 = wsumneg2 -sdif**2
                     wsectneg = wsumneg/ncalls
                     wvarneg2 = (wsumneg2/ncalls - wsectneg**2)/ncalls
                  endif
                  sumsez(nphot+1) = sumsez(nphot+1) - sdif
                  sumsez2(nphot+1) = sumsez2(nphot+1) - sdif**2
                  sdif   = 0.d0
               endif
               idmpbck = idumperrors
               idumperrors = 1
               call dump_for_debug(1,p3,p4,qph)
               idumperrors = idmpbck
            endif
            varbefore = varw 

* Hit or miss
            call wraprng(rnd,1)
            csi = 1.d0*rnd(1)*tetto
            ratio = sdif/tetto
            if (ratio.gt.ratio_max) then 
               ratio_max = ratio
            endif
            if (sdif.gt.tetto) then
               i_over = i_over + 1
               wsumov  = wsumov  + sdif - tetto
               wsumov2 = wsumov2 + (sdif - tetto)**2
               wsectov = wsumov/ncalls
               wvarov2 = (wsumov2/ncalls - wsectov**2)/ncalls
               if (i_over.eq.1) then
                  call messages(4)
               endif
               call dump_for_debug(2,p3,p4,qph)
            endif
            iunww = 1
c--- unweight
c$$$            i_reg = 0
            abssdif = abs(sdif)
            if (csi.lt.abssdif) then
               if (sdif.lt.0.d0) then
                  iunww = -1
                  nevunwneg = nevunwneg + 1
               endif
               if (storing.eq.'y'.and.idumpweighted.eq.0)
     .          call storehoraceevent(sdif,p3,p4,qph)
    !.              call dumpevent(nweighted,sdif,p3,p4,qph)
               i_reg = 0
            endif
c-- end unweight
         enddo
         iev   = iev+1.d0*iunww
         if (iev.eq.nmax.and.idumpweighted.eq.0) istop = 1
* The output is written in the file 'outfile'
         if (iwrite.eq.1.or.istop.eq.1) then
            iwrite = 0
            ahit  = 1.d0*iev
            ashut = 1.d0*ncalls
            eff   = ahit/ashut
            xsect = tetto*eff
            var   = tetto*sqrt(abs(eff*(1.d0-eff))/ashut)
            open(11,file=noefile,status='unknown')
            read(11,*)nmax
            close(11)
            call print_output(iev,xsect,var,wsect,wvar2,eff,i_over,nneg,
     .           wsectneg,wvarneg2,wsectov,wvarov2,ncut)
           call writedistributions

        endif
* loop over events is stopped
      enddo

      if (storing.eq.'y') call finalizestorage
	
	
	
c      call printdist	
	
      print*,'GENERATION FINISHED. ALL FILES SAVED IN ',path
      stop
      end
