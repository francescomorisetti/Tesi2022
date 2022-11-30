      function elmat2formz(p1,p2,p3,p4,k,Qq,Qf)
      implicit none
      common/fordebugging/idebugging
      integer i,idebugging
      double precision dot,elmat2formz,q2
      complex*16 elmat2,ecomplex_,e_
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),k(0:3),
     > ptmp(0:3), q(0:3)
      double precision kp1,kp2,kp3,kp4,kp12,kp22,kp32,kp42
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision p1p22,p1p32,p1p42,p2p32,p2p42,p3p42
      double precision p1p1,p2p2,p3p3,p4p4
      double precision mf,mf2,mf4,mq,mq2,mq4,s,mzm2,mzm4
      double precision Qq,Qf,Qq2,Qf2,Qq4,Qf4
      double precision kp1m1,kp2m1,kp3m1,kp4m1,kp1m2,kp2m2,kp3m2,kp4m2
      complex*16 im
      double precision charges(4),amasses(4),xx,ombet
      complex*16 eps1234
      double precision gw,gz,s12mmz2,s34mmz2,dg12,dg122,dg34,dg342
      double precision dz12m2,dz34m2
      double precision e2,e4,gs2sctw2,gs2sctw4,g2,gzmz,gzmz2
      double precision gaf,gvf,gaq,gvq,gaf2,gvf2,gaq2,gvq2,gfp
      double precision alpha,alpha_gf,s2tw,gf,pi,convfac
      double precision m1,m2,ch1,ch2,chfs,mfs
      double precision meorig,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt
      integer lepton,imirror
      common/t_widths/gw,gz
      double precision dum1(-5:5),dum2(-5:5),dum3(-5:5)
      common/masseshor/meorig,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     dum1,dum2,dum3
      common/processhor/m1,m2,ch1,ch2,chfs,mfs,lepton,imirror
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/chargesmasses_cm/charges,amasses
      common/imaginaryzrad/im
      data im /(0.d0,1.d0)/

      e2 = alpha*4.d0*pi
      e4 = e2*e2
      g2 = e2/s2tw
c      gfp = sqrt(2.d0)*e2/s2tw/8.d0/mw/mw
c      g2 = gfp/sqrt(2.d0)*mw*mw*8.d0
      gs2sctw2 = g2/4.d0/(1.d0-s2tw)
      gs2sctw4 = gs2sctw2*gs2sctw2

      gaf  = -0.5d0 ! electron weak isospin!
      gvf  =  gaf - 2.d0*Qf*s2tw
      gaf2 =  gaf*gaf
      gvf2 =  gvf*gvf
      if (Qq.gt.0) then
         gaq = 0.5d0
         gvq = gaq - 2.d0*Qq*s2tw
      else
         gaq = -0.5d0
         gvq =  gaq - 2.d0*Qq*s2tw
      endif
      gaq2 = gaq*gaq
      gvq2 = gvq*gvq

      Qf2 =Qf*Qf
      Qf4 =Qf2*Qf2
      Qq2 =Qq*Qq
      Qq4 =Qq2*Qq2

      mf  = amasses(3)
      mf2 = mf*mf
      mf4 = mf2*mf2
      mq  = amasses(1)
      mq2 = mq*mq

      kp1 = dot(k,p1)
      kp2 = dot(k,p2)
      kp3 = dot(k,p3)
      kp4 = dot(k,p4)

c      print*,'sqz.f line 69: curare kp3 e kp4!'
cc ! trying to avoid NaNs at extremely high energies      
      if (kp3.eq.0.d0) then
         xx    = -mf2/p3(0)/p3(0)
         ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .        - 0.078125d0 * xx*xx*xx)
         kp3 = p3(0)*k(0)*ombet
      endif
      if (kp4.eq.0.d0) then
         xx    = -mf2/p4(0)/p4(0)
         ombet = - 0.5d0*xx*(1.d0 - xx * 0.250d0 + xx*xx * 0.125d0
     .        - 0.078125d0 * xx*xx*xx)
         kp4 = p4(0)*k(0)*ombet
      endif
cc
      
      p1p1 = dot(p1,p1)
      p2p2 = dot(p2,p2)
      p3p3 = dot(p3,p3)
      p4p4 = dot(p4,p4)

      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4)
      p3p4 = dot(p3,p4)

      kp12 = kp1*kp1
      kp22 = kp2*kp2
      kp32 = kp3*kp3
      kp42 = kp4*kp4
      
      p1p22 = p1p2*p1p2
      p1p32 = p1p3*p1p3
      p1p42 = p1p4*p1p4
      p2p32 = p2p3*p2p3
      p2p42 = p2p4*p2p4
      p3p42 = p3p4*p3p4

      mq4 = mq2*mq2

      mzm2 = 1.d0/mz/mz
      mzm4 = mzm2*mzm2
      gzmz = gz*mz
      gzmz2= gzmz*gzmz

      do i = 0,3
         ptmp(i) = p1(i) + p2(i)
      enddo      
      q2      = dot(ptmp,ptmp)
      dg12    = 1.d0/q2
      dg122   = dg12*dg12
      s12mmz2 = q2 - mz*mz
      dz12m2  = 1.d0/((q2-mz**2)**2 + gzmz2)

      do i = 0,3
         ptmp(i) = p4(i) + p3(i)
      enddo
      q2      = dot(ptmp,ptmp)
      dg34    = 1.d0/q2
      dg342   = dg34*dg34
      s34mmz2 = q2 - mz*mz
      dz34m2  = 1.d0/((q2-mz**2)**2 + gzmz2)

      kp1m1 = 1.d0/kp1
      kp2m1 = 1.d0/kp2
      kp3m1 = 1.d0/kp3
      kp4m1 = 1.d0/kp4

      kp1m2 = kp1m1*kp1m1
      kp2m2 = kp2m1*kp2m1
      kp3m2 = kp3m1*kp3m1
      kp4m2 = kp4m1*kp4m1

****      eps1234 = e_(p1,p2,p3,p4)
      eps1234 = (0.d0,0.d0)
** Commenti: 
* con eps1234 acceso ho verificato che la parte immaginaria 
* dell'el. di m. e' O(10^-15), come deve essere. Se metto eps1234 = 0, 
* la parte imm. e' zero secco, ma cambia un po' anche la parte reale.
* La gauge invarianza e' checkata perfettamente, anche con km,kn a zero
** PER QUESTI TEST USO km,kn = 0 e eps1234 = 0.d0!!!

      elmat2 = (0.d0,0.d0)
      include 'form/zrad.f'
      elmat2formz = elmat2
      return
      end

***********************************
      function elmat2bornz(p1,p2,p3,p4,iba)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),ptmp(0:3)
      integer i,iba
      double precision dot,elmat2bornz,shardwrapper
      external dot,shardwrapper
      complex*16 elmat2
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision p1p22,p1p32,p1p42,p2p32,p2p42,p3p42
      double precision p1p1,p2p2,p3p3,p4p4
      double precision meorig,tmp1,tmp2,tmp3,tmp4
      complex*16 im
      double precision mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt
      double precision mf,mf2,mf4,mq,mq2,mq4,Qf,Qf2,Qf4,Qq,Qq2,Qq4
      double precision gw,gz,e2,e4,Qu,Qd,Qe,Qu2,Qd2,Qe2,ovall,xsw2l,
     .     xsw2u,xsw2d,s2twlepexp
      double precision gs2sctw2,gs2sctw4,mzm2,mzm4,deng,deng2,s,denzmod2
      complex*16 denz,denzd
      double precision gaf,gvf,gaq,gvq,gaf2,gvf2,gaq2,gvq2,g2,gfp
      double precision alpha,alpha_gf,s2tw,gf,pi,convfac,rhofi,dtmp
      double precision charges(4),amasses(4)
      character*5 zIS
      common/selectvirtualz/zIS
      common/t_widths/gw,gz
      double precision dum1(-5:5),dum2(-5:5),dum3(-5:5)
      common/masseshor/meorig,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     dum1,dum2,dum3
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/sin2twefflept/s2twlepexp
      common/chargesmasses_cm/charges,amasses
      common/imaginary/im
      data im /(0.d0,1.d0)/

      integer ischeme
      double precision dalpha,dalphaan,drho,drhoir1,drhoir12,rhofiu,
     .     rhofid,minnie, zeta, xsw2, cw2,dkl,dku,dkd
      common/feyhorflags/dalpha,dalphaan,drho,drhoir1,drhoir12,rhofiu,
     .     rhofid,dkl,dku,dkd,ischeme


      xsw2l = s2tw
      xsw2u = s2tw
      xsw2d = s2tw
      xsw2  = s2tw

      ovall = 32.d0

      e2 = alpha*4.d0*pi
      e4 = e2*e2
      g2 = e2/s2tw

c      gfp = sqrt(2.d0)*e2/s2tw/8.d0/mw/mw
c      g2 = gfp/sqrt(2.d0)*mw*mw*8.d0

      gs2sctw2 = g2/4.d0/(1.d0-s2tw)
      gs2sctw4 = gs2sctw2*gs2sctw2

      
      if (ischeme.gt.1.and.iba.eq.1) then
c versione "linearizzata"
c         minnie = 1d0+dalpha
c         e2 = alpha*4.d0*pi*minnie
c         e4 = (alpha*4.d0*pi)**2 * (1.d0 + 2.d0*dalpha)
         minnie = 1d0/(1d0-dalpha)
         e2 = alpha*4.d0*pi*minnie
         e4 = e2*e2            
         
         rhofi = rhofid
         if (zIS.eq.'uubar') rhofi = rhofiu

         dtmp = rhofi-1.d0

         gs2sctw2 = sqrt(2d0)*gf*mz*mz * rhofi/(1.d0-drhoir12)
         gs2sctw4 = gs2sctw2*gs2sctw2

c versione "linearizzata"         
c         gs2sctw2=sqrt(2d0)*gf*mz*mz * (1.d0+dtmp+drhoir1)
c         gs2sctw4=(sqrt(2d0)*gf*mz*mz)**2*(1.d0+2.d0*dtmp+2.d0*drhoir1)

         if (ischeme.eq.2) then
c            xsw2l = s2twlepexp
c            xsw2u = s2twlepexp - 0.0001d0
c            xsw2d = s2twlepexp - 0.0002d0
            xsw2l = s2tw*(1.d0+dkl)
            xsw2u = s2tw*(1.d0+dku)
            xsw2d = s2tw*(1.d0+dkd)
         endif
      endif

      Qf  = charges(3)
      Qq  = charges(1)
      Qf2 = Qf*Qf
      Qf4 = Qf2*Qf2
      Qq2 = Qq*Qq
      Qq4 = Qq2*Qq2

      gaf  = -0.5d0 ! electron weak isospin!
      gvf  =  gaf - 2.d0*charges(3)*xsw2l
      gaf2 =  gaf*gaf
      gvf2 =  gvf*gvf
      if (Qq.gt.0) then
         gaq = 0.5d0
         gvq = gaq - 2.d0*charges(1)*xsw2u
      else
         gaq = -0.5d0
         gvq =  gaq - 2.d0*charges(1)*xsw2d
      endif
      gaq2 = gaq*gaq
      gvq2 = gvq*gvq

      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4)
      p3p4 = dot(p3,p4)
      
      p1p22 = p1p2*p1p2
      p1p32 = p1p3*p1p3
      p1p42 = p1p4*p1p4
      p2p32 = p2p3*p2p3
      p2p42 = p2p4*p2p4
      p3p42 = p3p4*p3p4

      mf  = amasses(3)
      mq  = amasses(1)
      mf2 = mf*mf
      mq2 = mq*mq
c      mf2 = 0.d0
c      mq2 = 0.d0
      mf4 = mf2*mf2
      mq4 = mq2*mq2

      mzm2 = 1.d0/mz/mz
      mzm4 = mzm2*mzm2

      do i = 0,3
         ptmp(i) = p1(i) + p2(i)
      enddo      
      s = dot(ptmp,ptmp)

      deng  = 1.d0/s
      deng2 = deng*deng

      denzmod2 = 1.d0/((s-mz**2)**2 + gz**2*mz**2)
      denz     = 1.d0/(s-mz**2 + im*gz*mz)
      denzd    = 1.d0/(s-mz**2 - im*gz*mz)

      elmat2 = (0.d0,0.d0)
      include 'form/bornz.f'
      elmat2bornz = 1.d0*elmat2
      return
      end
c===========
      function elmat2borngg(p1,p2,p3,p4)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),ptmp(0:3)
      integer i
      double precision dot,elmat2borngg
      external dot
      complex*16 elmat2
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision p1p22,p1p32,p1p42,p2p32,p2p42,p3p42
      double precision p1p1,p2p2,p3p3,p4p4
      double precision p1p3m2,p2p3m2,p1p3m1,p2p3m1
      double precision meorig,tmp1,tmp2,tmp3,tmp4,me
      double precision mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt
      double precision mf,mf2,mf4,mq,mq2,mq4,Qf,Qf2,Qf4,Qq,Qq2,Qq4
      double precision gw,gz,e2,e4,Qu,Qd,Qe,Qu2,Qd2,Qe2,ovall
      double precision gs2sctw2,gs2sctw4,mzm2,mzm4,deng,deng2,s,denzmod2
      double precision alpha,alpha_gf,s2tw,gf,pi,convfac
      double precision charges(4),amasses(4)
      double precision dum1(-5:5),dum2(-5:5),dum3(-5:5)
      common/masseshor/meorig,mmu,mw,mz,mu,md,mc,ms,mtau,mb,mt,
     .     dum1,dum2,dum3
      common/const/alpha,alpha_gf,s2tw,gf,pi,convfac
      common/chargesmasses_cm/charges,amasses

      e2 = alpha*4.d0*pi
      e4 = e2*e2

      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4)
      p3p4 = dot(p3,p4)

      p1p22 = p1p2*p1p2
      p1p32 = p1p3*p1p3
      p1p42 = p1p4*p1p4
      p2p32 = p2p3*p2p3
      p2p42 = p2p4*p2p4
      p3p42 = p3p4*p3p4


      p1p3m2 = 1.d0/p1p32
      p2p3m2 = 1.d0/p2p32
      p1p3m1 = 1.d0/p1p3
      p2p3m1 = 1.d0/p2p3

      mf  = amasses(3)
      mf2 = mf*mf
      mf4 = mf2*mf2

      elmat2 = (0.d0,0.d0)
      include 'form/ggborn.f'
      elmat2borngg = 1.d0*elmat2

      return
      end
***
      function shardwrapper(p1,p2,p3,p4)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),ptmp(0:3)
      common/partons/ipart1,ipart2
      do i = 0,3
         ptmp(i) = p1(i) + p2(i)
      enddo      
      s = dot(ptmp,ptmp)
      c = tridot(p1,p3)/sqrt(tridot(p3,p3))/sqrt(tridot(p1,p1))
      shardwrapper = shard(s,c,ipart1,ipart2)
      return
      end
