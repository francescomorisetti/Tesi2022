      implicit double precision (a-h,o-z)      
      parameter (ifl = 5) ! the total number of flavours playing in our
                          ! processes...
      double precision me,mmu,mw,mz,mu,md,mc,ms,mtau,m1,m2,mq(-ifl:ifl),
     >       chq(-ifl:ifl),iplus,mfs,itre(-ifl:ifl),mb,mt,mh
      double precision ckm(3,3),pdfpdf(-ifl:ifl,-ifl:ifl),
     >       myckm(-ifl:ifl,-ifl:ifl)!,sezhl(-ifl:ifl,-ifl:ifl)
      double precision pdf1(-ifl:ifl),pdf2(-ifl:ifl),p1(0:3),
     >    p2(0:3),p3(0:3),p4(0:3),qph(40,0:3),qph1(0:3),qph2(0:3),
     >    ptmp(0:3)
      integer hadr1,hadr2,charge_factor(4)
      character*5 ptqed,qedorder,horaceversion
      double precision nmax,iev,mass2
      character*1   boson
      character*1   storing
      character*100 outfile,readfile,noefile
      character*30  subtractionfile
      character*100 binfile,cutsfile
      character*150 storfile,conffile,command
      character*100 path
      character*72  line
      character*50  commentonpdfs
      double precision horaces,horacet,horaceu
      common/kinvars/horaces,horacet,horaceu
      common/masseshor/me,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,mq,chq,itre
      common/higgsmass/mh
      common/t_widths/gw,gz
      common/p_widths/gwe,gwm,gwtau,gze,gzm,gztau
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/sin2twefflept/s2twlepexp
      common/ckm_matrix/ckm,myckm
      common/is/hadr1,hadr2,i_pdf,iquark1,iquark2  ! initial state:
                                                   ! p - p, p - p_bar or
                                                   ! iquark1 - iquark2 (for
                                                   ! parton level x-section)
      common/photoninducedprocs/iphinduced
c      iphinduced = 1 ! 0 --> NO PHOTON INDUCED
c      iphinduced = 0 !--> NO PHOTON INDUCED

      common/pdf/pdf1,pdf2,pdfpdf
      common/pdfcomment/commentonpdfs
      common/processhor/m1,m2,ch1,ch2,chfs,mfs,lepton,imirror ! Stores 
                                                   ! the masses and 
						   ! charges of particles
                                                   ! involved in the process
      common/ichfinalstate/ichfs
      common/partons/ipart1,ipart2
      common/partonspdgmcscheme/ip1mc,ip2mc,ip3mc,ip4mc,ip5mc
      common/ps_param/eps,iplus
      common/exp_cuts/etamax,ptmin,cmin,cmax,ptminmiss,
     .     tmmin,mass2,ismear,irec,icutmt,ischemealpha
      common/pdf_parameters/qpdfmin,qpdfmax,xpdfmax,xpdfmin,
     >     rescalepdfscale,spdfcomm,iwhich,iscalechoice
      common/ch_factor/charge_factor
      common/QEDparam/ptqed,qedorder
      common/gen_event/ebeam,ebeam1,ebeam2,xmin,xmax,shat,x1pdf,x2pdf
*******************
! commons for remaining input variables
      common/vectorboson/boson
      common/int_input/nmax,iseed,iwriteoutput,idumpweighted
      common/real_input/bmin,bmax
      common/char_input/path,storing,outfile,noefile,
     >                  storfile,conffile,cutsfile,subtractionfile
      common/tunehorm/tune_hm,ratio_max
      common/new_old_kine/new_kine
      common/switcharun/ialphachoice
      common/iplsubtraction/iplsub
****
! common useful for massive quarks....
      common/quarksmassivevars/E1cm,E2cm,E1l,E2l,pz1,pz2,pz1l,pz2l
      common/x_pdf_scheme/i_scheme
**** 
! common for exact EWK...
      common/iexactewk/iewk
      common/sezionidiff/sezioni(20),sezioni2(20),sumsez(20),
     >     sumsez2(20),seztot(20),sezerr(20),sezij(-ifl:ifl,-ifl:ifl),
     >     sezij2(-ifl:ifl,-ifl:ifl),
     >     nphotmax,icounterpartial(20)
! HORACE VERSION
      common/horaceversion/horaceversion

! PDF REPLICA
      character*3 cpdfreplica
      common/pdfreplica/ipdfreplica,ilrep,cpdfreplica
