*******************************************
      subroutine initstorage
      include 'shared.inc'
      common/storageunit/nunitstor

      nutitstor = 66
      
! set ififo = 1 to get a bzip2-compressed event file on the fly
      ififo = 0
      
      if (ififo.eq.1) then  
         ilens = lnblnk(storfile)-5
         call
     >    EXECUTE_COMMAND_LINE("mkfifo "//storfile(1:ilens)//".fifo")
         
         call EXECUTE_COMMAND_LINE("bzip2 -9 < "//storfile(1:ilens)//
     >        ".fifo > "//storfile(1:lnblnk(storfile))//".bz2",
     >        wait=.false.)
         
         open
     >    (nunitstor,file=storfile(1:ilens)//".fifo",status='unknown')
         
      else
         open(nunitstor,file=storfile,status='unknown')
      endif

c     write(nunitstor,*)'Events produced by Horace '//horaceversion
c     write(nunitstor,*)'-------------------------------'
c     if (ichfs.eq.-1) write(nunitstor,*)'W- production'
c     if (ichfs.eq. 1) write(nunitstor,*)'W+ production'
c     if (ichfs.eq. 2) write(nunitstor,*)'W+ and W- production'
c     if (ichfs.eq. 0) write(nunitstor,*)'Z/gamma production'
c     write(nunitstor,*)'CoM energy',2.d0*ebeam
c     write(nunitstor,*)'Initial state hadrons',hadr1,hadr2
c     if (idumpweighed.eq.1) then
c     write(nunitstor,*)'weighted events'
c     else
c     write(nunitstor,*)'unweighted events'
c     endif
c     write(nunitstor,*)'-------------------------------'
      
      return
      end
*******************************************
      subroutine finalizestorage
      include 'shared.inc'
      common/storageunit/nunitstor

      close(nunitstor)
      
      end
*******************************************
      subroutine writeconffile(xsecup,xerrup,xmaxup)
      include 'shared.inc'
      common/nentries/nev
      nunit = 44
      open(nunit,file=conffile,status='unknown')
      ivr = lnblnk(horaceversion)
      write(nunit,*)'Horace '//horaceversion(1:ivr)
     $     //' LHA configuration'
      write(nunit,*)'----------------------------'
      write(nunit,*)'file where the events are stored'
      ilp = lnblnk(path)
      write(nunit,*)storfile(ilp+1:85)
      if (ichfs.eq.-1) write(nunit,*)'W- production'
      if (ichfs.eq. 1) write(nunit,*)'W+ production'
      if (ichfs.eq. 2) write(nunit,*)'W+ and W- production'
      if (ichfs.eq. 0) write(nunit,*)'Z/gamma production'
      write(nunit,*)commentonpdfs
      write(nunit,*)'CoM energy'
      write(nunit,*)2.d0*ebeam
      write(nunit,*)'Initial state hadrons'
      if (i_pdf.eq.1) then
         write(nunit,*)hadr1
         write(nunit,*)hadr2
      else
         write(nunit,*)iquark1
         write(nunit,*)iquark2
      endif
      if (idumpweighted.eq.1) then
         write(nunit,*)'weighted events'
         write(nunit,*)'idwtup'
         write(nunit,*)'-1'
      else
         write(nunit,*)'unweighted events'
         write(nunit,*)'idwtup'
         write(nunit,*)'-3'
      endif
      write(nunit,*)'xsecup'
      write(nunit,*)xsecup
      write(nunit,*)'xerrup'
      write(nunit,*)xerrup
      write(nunit,*)'xmaxup'
      write(nunit,*)xmaxup
      write(nunit,*)'1/alpha_em'
      write(nunit,*)1.d0/alpha
      write(nunit,*)'number of entries in the file (at least)'
      write(nunit,*)nev
      close(nunit)
      return
      end
c-------------------------------------------------------------------
      subroutine storehoraceevent(sez,p3,p4,qph)
c modified to dump in the same format as ALPGEN
c     dump event details to a file, for future reading by herwig
c-------------------------------------------------------------------
      include 'shared.inc'
      integer maxdec
      parameter (maxdec=40)
      dimension pin1(0:3),pin2(0:3)
      integer Nunit,npart,iproc,ifla(maxdec),icolup(2,maxdec)
      real*8 qsq,p(5,maxdec),wgt
      real Sq,Sp(3,maxdec),Sm(maxdec),Swgt
      integer nev,i,j
      common/nentries/nev
      common/storageunit/nunitstor
      common/evdumpHORACEcmn/ifirstdump
      data ifirstdump /0/
      data nev/0/
c
      if (ifirstdump.eq.0) then
         ifirstdump = 1
      endif

      nev = nev+1
      qsq = spdfcomm**2
      if (idumpweighted.eq.1) then
         wgt = sez
      else
         wgt = sez/abs(sez)
      endif
      nphot = 0
      do k = 1,maxdec
         if (qph(k,0).gt.0.d0) nphot = nphot + 1
      enddo
      npart=nphot + 4

      pin1(0) =  x1pdf * ebeam
      pin1(1) =  0.d0
      pin1(2) =  0.d0
      pin1(3) =  x1pdf * ebeam
      pin2(0) =  x2pdf * ebeam
      pin2(1) =  0.d0
      pin2(2) =  0.d0
      pin2(3) = -x2pdf * ebeam

      p(1,1)  = pin1(1)
      p(2,1)  = pin1(2)
      p(3,1)  = pin1(3)
      p(4,1)  = pin1(0)
      p(5,1)  = 0.d0

      p(1,2)  = pin2(1)
      p(2,2)  = pin2(2)
      p(3,2)  = pin2(3)
      p(4,2)  = pin2(0)
      p(5,2)  = 0.d0

      p(1,3)  = p3(1)
      p(2,3)  = p3(2)
      p(3,3)  = p3(3)
      p(4,3)  = p3(0)
      p(5,3)  = mfs

      p(1,4)  = p4(1)
      p(2,4)  = p4(2)
      p(3,4)  = p4(3)
      p(4,4)  = p4(0)
      p(5,4)  = 0.d0
      if (boson.eq.'Z') p(5,4) = mfs

      if (nphot.gt.0) then
         do k = 1,nphot
            p(1,k+4)  = qph(k,1)
            p(2,k+4)  = qph(k,2)
            p(3,k+4)  = qph(k,3)
            p(4,k+4)  = qph(k,0)
            p(5,k+4)  = 0.d0
            ifla(k+4) = ip5mc
            icolup(1,k+4) = 0
            icolup(2,k+4) = 0
            if (ip5mc.ne.22) then
               if (ip5mc.gt.0) then
                  icolup(1,k+4) = 501
                  icolup(2,k+4) = 0
               else
                  icolup(1,k+4) = 0
                  icolup(2,k+4) = 501
               endif
            endif
         enddo
      endif

      iproc   = 1
      ifla(1) = ip1mc
      ifla(2) = ip2mc
      ifla(3) = ip3mc
      ifla(4) = ip4mc

      icolup(1,1) = 0
      icolup(2,1) = 0
      icolup(1,2) = 0
      icolup(2,2) = 0
      if (ip1mc.ne.22) then
         if (ip1mc.gt.0) then 
            icolup(1,1) = 501
            icolup(2,1) = 0
         else
            icolup(1,1) = 0
            icolup(2,1) = 501
         endif
      endif
      if (ip2mc.ne.22) then
         if (ip2mc.gt.0) then 
            icolup(1,2) = 501
            icolup(2,2) = 0
         else
            icolup(1,2) = 0
            icolup(2,2) = 501
         endif
      endif
      icolup(1,3) = 0
      icolup(2,3) = 0
      icolup(1,4) = 0
      icolup(2,4) = 0

      Sq   = real(sqrt(qsq)) ! PDF scale
      Swgt = real(wgt)       ! weight, +-1

      do i=1,npart
         do j=1,3
            Sp(j,i) = real(p(j,i))  ! 1,2,3 --> x,y,z
         enddo
         Sm(i) = real(p(5,i))       ! mass
      enddo

***      if (npart.eq.5.and.ip5mc.ne.22) then
c     Nevent, iproc, npart, wgt, Q scale 
      write(nunitstor,2) '>',nev,npart,Swgt,Sq  ! iproc, random.... (1)
 2    format((A),i8,1x,i4,1x,2(1x,e12.6))
c     flavour, colour and z-momentum of incoming partons
      write(nunitstor,8) ifla(1),icolup(1,1),icolup(2,1),Sp(3,1) ! ifl, ....
      write(nunitstor,8) ifla(2),icolup(1,2),icolup(2,2),Sp(3,2)
c     flavour, colour, 3-momentum and mass of outgoing partons
      do i=3,npart
         write(nunitstor,9) ifla(i),icolup(1,i),icolup(2,i),
     $        Sp(1,i),Sp(2,i),Sp(3,i),Sm(i)
      enddo
 8    format(i8,1x,2(i4,1x),f12.5)
 9    format(i8,1x,2(i4,1x),4(1x,f12.5))
***      endif
      return
      end
**************************************************************
      subroutine evdumpHORACEglasgow(sez,p3,p4,qph)
c used for the "Physics with early LHC data" workshop at Glasgow
c modified to dump in the same format as ALPGEN
c     dump event details to a file, for future reading by herwig
c-------------------------------------------------------------------
      include 'shared.inc'
      integer maxdec
      parameter (maxdec=40)
      dimension pin1(0:3),pin2(0:3)
      integer Nunit,npart,iproc,ifla(maxdec),icolup(2,maxdec)
      real*8 qsq,p(5,maxdec),wgt
      real Sq,Sp(3,maxdec),Sm(maxdec),Swgt
      integer nev,i,j
      common/evdumpHORACEglasgowcmn/ifirstdump
      data ifirstdump /0/
      data nev/0/
      save
c
      Nunit = 66
      if (ifirstdump.eq.0) then
         open (Nunit,file=storfile,status='unknown')
      endif
      nev = nev+1
      qsq = spdfcomm**2
      if (idumpweighted.eq.1) then
         wgt = sez
      else
         wgt = sez/abs(sez)
      endif
      nphot = 0
      do k = 1,maxdec
         if (qph(k,0).gt.0.d0) nphot = nphot + 1
      enddo
      npart=nphot + 4

      pin1(0) =  x1pdf * ebeam
      pin1(1) =  0.d0
      pin1(2) =  0.d0
      pin1(3) =  x1pdf * ebeam
      pin2(0) =  x2pdf * ebeam
      pin2(1) =  0.d0
      pin2(2) =  0.d0
      pin2(3) = -x2pdf * ebeam

      p(1,1)  = pin1(1)
      p(2,1)  = pin1(2)
      p(3,1)  = pin1(3)
      p(4,1)  = pin1(0)
      p(5,1)  = 0.d0

      p(1,2)  = pin2(1)
      p(2,2)  = pin2(2)
      p(3,2)  = pin2(3)
      p(4,2)  = pin2(0)
      p(5,2)  = 0.d0

      p(1,3)  = p3(1)
      p(2,3)  = p3(2)
      p(3,3)  = p3(3)
      p(4,3)  = p3(0)
      p(5,3)  = mfs

      p(1,4)  = p4(1)
      p(2,4)  = p4(2)
      p(3,4)  = p4(3)
      p(4,4)  = p4(0)
      p(5,4)  = 0.d0
      if (boson.eq.'Z') p(5,4) = mfs

      if (nphot.gt.0) then
         do k = 1,nphot
            p(1,k+4)  = qph(k,1)
            p(2,k+4)  = qph(k,2)
            p(3,k+4)  = qph(k,3)
            p(4,k+4)  = qph(k,0)
            p(5,k+4)  = 0.d0
            ifla(k+4) = ip5mc
            icolup(1,k+4) = 0
            icolup(2,k+4) = 0
            if (ip5mc.ne.22) then
               if (ip5mc.gt.0) then
                  icolup(1,k+4) = 501
                  icolup(2,k+4) = 0
               else
                  icolup(1,k+4) = 0
                  icolup(2,k+4) = 501
               endif
            endif
         enddo
      endif

      iproc   = 1
      ifla(1) = ip1mc
      ifla(2) = ip2mc
      ifla(3) = ip3mc
      ifla(4) = ip4mc

      icolup(1,1) = 0
      icolup(2,1) = 0
      icolup(1,2) = 0
      icolup(2,2) = 0
      if (ip1mc.ne.22) then
         if (ip1mc.gt.0) then 
            icolup(1,1) = 501
            icolup(2,1) = 0
         else
            icolup(1,1) = 0
            icolup(2,1) = 501
         endif
      endif
      if (ip2mc.ne.22) then
         if (ip2mc.gt.0) then 
            icolup(1,2) = 501
            icolup(2,2) = 0
         else
            icolup(1,2) = 0
            icolup(2,2) = 501
         endif
      endif
      icolup(1,3) = 0
      icolup(2,3) = 0
      icolup(1,4) = 0
      icolup(2,4) = 0

      Sq   = real(sqrt(qsq)) ! PDF scale
      Swgt = real(wgt)       ! weight, +-1

      do i=1,npart
         do j=1,3
            Sp(j,i) = real(p(j,i))  ! 1,2,3 --> x,y,z
         enddo
         Sm(i) = real(p(5,i))       ! mass
      enddo

***      if (npart.eq.5.and.ip5mc.ne.22) then
c     Nevent, iproc, npart, wgt, Q scale 
      iii = 1
      write(Nunit,2) nev,iii,npart,Swgt,Sq  ! iproc, random.... (1)
 2    format(i8,1x,i4,i3,2(1x,e12.6))
c     flavour, colour and z-momentum of incoming partons
      write(Nunit,8) ifla(1),icolup(1,1),icolup(2,1),Sp(3,1) ! ifl, ....
      write(Nunit,8) ifla(2),icolup(1,2),icolup(2,2),Sp(3,2)
c     flavour, colour, 3-momentum and mass of outgoing partons
      do i=3,npart
         write(Nunit,9) ifla(i),icolup(1,i),icolup(2,i),Sp(1,i),Sp(2,i),
     $        Sp(3,i),Sm(i)
      enddo
 8    format(i8,1x,2(i4,1x),f10.3)
 9    format(i8,1x,2(i4,1x),4(1x,f10.3))
***      endif
      return
      end
c-------------------------------------------------------------------
      subroutine dumpevent(nc,sez,p3,p4,qph)
      include 'shared.inc'
      integer*8 nc
      common/dumpeventcommon/ifirstdump,icount
      data ifirstdump,icount /0,0/
      if (ifirstdump.eq.0) then

         print*,'This event format is too naif!!'
         stop

         open (66,file=storfile,status='unknown')
         if (idumpweighted.eq.1) then
            write(66,*)'Weighted events produced by HORACE ',
     .           horaceversion
         else
            write(66,*)'Unweighted events produced by HORACE ',
     .           horaceversion
         endif
         write(66,*)nmax
         write(66,*)2.d0*ebeam
         write(66,*)hadr1
         write(66,*)hadr2
         write(66,*)alpha
         ifirstdump = 1
      endif
      icount = icount + 1
      write(66,*)icount,' -----------'
      write(66,*)nc
      nphot = 0
      do k = 1,40
         if (qph(k,0).gt.0.d0) nphot = nphot + 1
      enddo
      write(66,*)nphot
      write(66,*)x1pdf,x2pdf,spdfcomm
      write(66,*)ipart1,ipart2
      if (lepton.eq.1.and.chfs.gt.0.d0) then
         ipart3 = -11
         ipart4 = 12
      endif
      if (lepton.eq.1.and.chfs.lt.0.d0) then
         ipart3 =  11
         ipart4 = -12
      endif
      if (lepton.eq.2.and.chfs.gt.0.d0) then
         ipart3 = -13
         ipart4 = 14
      endif
      if (lepton.eq.2.and.chfs.lt.0.d0) then
         ipart3 =  13
         ipart4 = -14
      endif

      if (lepton.eq.3.and.chfs.gt.0.d0) then
         ipart3 = -15
         ipart4 = 16
      endif
      if (lepton.eq.3.and.chfs.lt.0.d0) then
         ipart3 =  15
         ipart4 = -16
      endif

      write(66,*)ipart3,ipart4
      if (idumpweighted.eq.1) then
         write(66,*)sez
      else
         write(66,*)sez/abs(sez)
      endif
      write(66,*)p3
      write(66,*)p4
      iprimi = 0
      if (iewk.eq.0) iprimi = 20
      if (nphot.gt.0) then
         do k = 1+iprimi,nphot+iprimi
            write(66,*)(qph(k,i),i=0,3)
         enddo
      endif
      return
      end
**********
      subroutine dump_for_debug(i,p3,p4,qph)
      include 'shared.inc'
      integer isvec(25)
      common/rlxisvec/isvec
      
      character*150 errorfile
      dimension p1cm(0:3),p2cm(0:3),p3cm(0:3),p4cm(0:3)
      dimension p1b(0:3),p2b(0:3)
      dimension qphcm(40,0:3)
      dimension qtmp(0:3)      
      common/momentainitial/p1cm,p2cm
      common/eventinCOM/p3cm,p4cm,qphcm
      common/reducedtoborn/p1b,p2b,iref
      common/momentainitialred/pin1b(0:3),pin2b(0:3),pin1r(0:3)
     >     ,pin2r(0:3)
      common/radpattern/nph(4)
      common/svfordebug/svdebug,emtx,emtxsub,wnphot
      common/n_phot/nphot
      common/phindchannel/bbb,pw0,qmod,c1q,ich
      common/switchdumperrors/idumperrors
      common/fordumping/varw,varbefore,iev,idestroyed
      common/propagateurs/pr1,pr2,pr3,pr4,pr5,pr6,pr7,pr8,pr9,pr10
      common/propstest/testprop1,testprop2
      common/firstdump/ifirst
      data ifirst /1/

      if (idumperrors.eq.0) return

c      call printvector(i,p3cm,p4cm,qphcm)

c      do k = 1,100
c         if (path(k:k).ne.' ') lpath=k
c      enddo

      
c      errorfile=path(:lpath)//'errors.txt'

      iii = lnblnk(outfile)

      errorfile = storfile(1:iii-4)//'err'
      
      if (ifirst.eq.1) then
         open(99,file=errorfile,status='unknown')
      else
         open(99,file=errorfile,status='unknown',access='append')
      endif
      write(99,*)'====== DUMPING at point ',i
      if (idestroyed.eq.1) then
         write(99,*)'  this event was rejected!!'
      endif
      write(99,*)'propag ',pr1,pr2,pr3,pr4
      write(99,*)'test props. ',testprop1,testprop2
      write(99,*)'Var        ',varw
      write(99,*)'Var before ',varbefore
      write(99,*)'iev, nphot and rad. pattern ',iev,nphot,'   ',nph
      write(99,*)'partons... ',ipart1,ipart2
      write(99,*)'x1 & x2 ',x1pdf,x2pdf
      write(99,*)'sv & ME ',svdebug,emtx
      write(99,*)'plepton   = ',p3cm
      write(99,*)'pneutrino = ',p4cm
      write(99,*)'pleptonlab   = ',p3
      write(99,*)'pneutrinolab = ',p4
      do k = 0,3
         ptmp(k) = p3cm(k) + p4cm(k)
      enddo
      sqrtq2 = sqrt(dot(ptmp,ptmp))

      do k = 0,3
         ptmp(k) = p1cm(k) + p2cm(k)
      enddo
      write(99,*)'sqrt(s) = ',sqrt(dot(ptmp,ptmp))
      write(99,*)'sqrt_fs_cm = ', sqrtq2
      write(99,*)'phind channel', ich
      write(99,*)'bbb pw(0) qmod c1q',bbb,pw0,qmod,c1q
      write(99,*)'reduced to Born momenta'
      write(99,*)pin1b
      write(99,*)pin2b
      write(99,*)p1b
      write(99,*)p2b
      if (nphot.gt.0) then
         do j = 1,nphot
         write(99,*)'>> Photon n.',j
         do k = 0,3
            ptmp(k) = qphcm(j,k)            
            qtmp(k) = qph(j,k)
         enddo
         write(99,*)'cm momentum  ',ptmp
         write(99,*)'lab momentum ',qtmp
         write(99,*)'c_ph_up  ',tridot(p1cm,ptmp)/
     .        sqrt(tridot(p1cm,p1cm))/ptmp(0)
         write(99,*)'c_ph_up lab ',qtmp(3)/max(qtmp(0),1.d-12)
         write(99,*)'c_ph_lep',tridot(p3cm,ptmp)/
     .        sqrt(tridot(p3cm,p3cm))/ptmp(0)
         write(99,*)'c_ph_lep lab ',tridot(p3,qtmp)/
     .        sqrt(tridot(p3,p3))/max(qtmp(0),1.d-12)
         enddo
      endif
      write(99,*)'isvec for ranlux '
      write(99,*)isvec

      write(99,*)' '
      close(99)      
      ifirst=0
      return
      end
*************
      subroutine eventstorage(p3,p4,qph)
      include 'shared.inc'
      parameter (nsize_event = 16)
      real*4 event(nsize_event) 

      return
      end
      subroutine init_chtags(chtags,n)
      integer n
      character*4 chtags(n)
      chtags(1)  = 'e1'
      chtags(2)  = 'p1x'
      chtags(3)  = 'p1y'
      chtags(4)  = 'p1z'
      chtags(5)  = 'e2'
      chtags(6)  = 'p2x'
      chtags(7)  = 'p2y'
      chtags(8)  = 'p2z'
      chtags(9)  = 'q10'
      chtags(10) = 'q1x'
      chtags(11) = 'q1y'
      chtags(12) = 'q1z'
      chtags(13) = 'q20'
      chtags(14) = 'q2x'
      chtags(15) = 'q2y'
      chtags(16) = 'q2z'
      return
      end
