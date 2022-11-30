      elmat2 =
     &  + mwc2dm1*md2 * (  - 64.D0*p1p4*p2p3 - 64.D0*p1p3*p2p4 - 128.D0
     &    *p1p3*p1p4 + 64.D0*p1p2*p3p4 )
      elmat2 = elmat2 + mwc2dm1*mup2 * (  - 128.D0*p2p3*p2p4 - 64.D0*
     &    p1p4*p2p3 - 64.D0*p1p3*p2p4 + 64.D0*p1p2*p3p4 )
      elmat2 = elmat2 + mwc2dm1*mup2*md2 * ( 128.D0*p3p4 )
      elmat2 = elmat2 + mwc2m1*md2 * (  - 64.D0*p1p4*p2p3 - 64.D0*p1p3*
     &    p2p4 - 128.D0*p1p3*p1p4 + 64.D0*p1p2*p3p4 )
      elmat2 = elmat2 + mwc2m1*mup2 * (  - 128.D0*p2p3*p2p4 - 64.D0*
     &    p1p4*p2p3 - 64.D0*p1p3*p2p4 + 64.D0*p1p2*p3p4 )
      elmat2 = elmat2 + mwc2m1*mup2*md2 * ( 128.D0*p3p4 )
      elmat2 = elmat2 + mwc2m1*mwc2dm1*md4 * (  - 64.D0*p1p2*p3p4 )
      elmat2 = elmat2 + mwc2m1*mwc2dm1*md2 * ( 128.D0*p1p2*p2p3*p2p4 + 
     &    128.D0*p1p2*p1p4*p2p3 + 128.D0*p1p2*p1p3*p2p4 + 128.D0*p1p2*
     &    p1p3*p1p4 - 128.D0*p1p22*p3p4 )
      elmat2 = elmat2 + mwc2m1*mwc2dm1*mup4 * (  - 64.D0*p1p2*p3p4 )
      elmat2 = elmat2 + mwc2m1*mwc2dm1*mup4*md2 * (  - 128.D0*p3p4 )
      elmat2 = elmat2 + mwc2m1*mwc2dm1*mup2 * ( 128.D0*p1p2*p2p3*p2p4
     &     + 128.D0*p1p2*p1p4*p2p3 + 128.D0*p1p2*p1p3*p2p4 + 128.D0*
     &    p1p2*p1p3*p1p4 - 128.D0*p1p22*p3p4 )
      elmat2 = elmat2 + mwc2m1*mwc2dm1*mup2*md4 * (  - 128.D0*p3p4 )
      elmat2 = elmat2 + mwc2m1*mwc2dm1*mup2*md2 * ( 256.D0*p2p3*p2p4 + 
     &    256.D0*p1p4*p2p3 + 256.D0*p1p3*p2p4 + 256.D0*p1p3*p1p4 - 384.D
     &    0*p1p2*p3p4 )
      elmat2 = elmat2 + 256.D0*p1p4*p2p3
