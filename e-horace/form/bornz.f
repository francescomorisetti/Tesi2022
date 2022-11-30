      elmat2 =                                                          
     &  + ovall*deng2*Qf2*Qq2*e4 * ( 2.D0*mq2*mf2 + p3p4*mq2 + p1p4*    
     &    p2p3 + p1p3*p2p4 + p1p2*mf2 )                                 
      elmat2 = elmat2 + ovall*gs2sctw4*denzmod2*gvf2*gvq2 * ( 2.D0*mq2* 
     &    mf2 + p3p4*mq2 + p1p4*p2p3 + p1p3*p2p4 + p1p2*mf2 )           
      elmat2 = elmat2 + ovall*gs2sctw4*denzmod2*gaq2*gvf2 * ( 2.D0*mq4* 
     &    mf2*mzm2 - 2.D0*mq2*mf2 + 2.D0*p3p4*mq4*mzm2 - p3p4*mq2 - 2.D0
     &    *p2p3*p2p4*mq2*mzm2 + p1p4*p2p3 - 2.D0*p1p4*p2p3*mq2*mzm2 +   
     &    p1p3*p2p4 - 2.D0*p1p3*p2p4*mq2*mzm2 - 2.D0*p1p3*p1p4*mq2*mzm2 
     &     + p1p2*mf2 + 2.D0*p1p2*mq2*mf2*mzm2 + 2.D0*p1p2*p3p4*mq2*    
     &    mzm2 )                                                        
      elmat2 = elmat2 + ovall*gs2sctw4*denzmod2*gaf2*gvq2 * ( 2.D0*mq2* 
     &    mf4*mzm2 - 2.D0*mq2*mf2 + p3p4*mq2 + 2.D0*p3p4*mq2*mf2*mzm2   
     &     - 2.D0*p1p4*p2p4*mf2*mzm2 + p1p4*p2p3 - 2.D0*p1p4*p2p3*mf2*  
     &    mzm2 + p1p3*p2p4 - 2.D0*p1p3*p2p4*mf2*mzm2 - 2.D0*p1p3*p2p3*  
     &    mf2*mzm2 + 2.D0*p1p2*mf4*mzm2 - p1p2*mf2 + 2.D0*p1p2*p3p4*mf2 
     &    *mzm2 )                                                       
      elmat2 = elmat2 + ovall*gs2sctw4*denzmod2*gaf2*gaq2 * (  - 2.D0*  
     &    mq4*mf2*mzm2 - 2.D0*mq2*mf4*mzm2 + 2.D0*mq2*mf2 + 2.D0*p3p4*  
     &    mq4*mzm2 - p3p4*mq2 - 2.D0*p3p4*mq2*mf2*mzm2 - 2.D0*p2p3*p2p4 
     &    *mq2*mzm2 + 4.D0*p2p3*p2p4*mq2*mf2*mzm4 - 2.D0*p1p4*p2p4*mf2* 
     &    mzm2 + 4.D0*p1p4*p2p4*mq2*mf2*mzm4 + p1p4*p2p3 - 2.D0*p1p4*   
     &    p2p3*mf2*mzm2 - 2.D0*p1p4*p2p3*mq2*mzm2 + 4.D0*p1p4*p2p3*mq2* 
     &    mf2*mzm4 + p1p3*p2p4 - 2.D0*p1p3*p2p4*mf2*mzm2 - 2.D0*p1p3*   
     &    p2p4*mq2*mzm2 + 4.D0*p1p3*p2p4*mq2*mf2*mzm4 - 2.D0*p1p3*p2p3* 
     &    mf2*mzm2 + 4.D0*p1p3*p2p3*mq2*mf2*mzm4 - 2.D0*p1p3*p1p4*mq2*  
     &    mzm2 + 4.D0*p1p3*p1p4*mq2*mf2*mzm4 + 2.D0*p1p2*mf4*mzm2 -     
     &    p1p2*mf2 - 2.D0*p1p2*mq2*mf2*mzm2 + 2.D0*p1p2*p3p4*mf2*mzm2   
     &     + 2.D0*p1p2*p3p4*mq2*mzm2 + 2.D0*p2p42*mq2*mf2*mzm4 + 2.D0*  
     &    p2p32*mq2*mf2*mzm4 + 2.D0*p1p42*mq2*mf2*mzm4 + 2.D0*p1p32*mq2 
     &    *mf2*mzm4 )                                                   
      elmat2 = elmat2 + ovall*gaf*gaq*gvf*gvq*gs2sctw4*denzmod2 * ( 4.D0
     &    *p1p4*p2p3 - 4.D0*p1p3*p2p4 )                                 
      elmat2 = elmat2 + e2*ovall*gvf*gvq*gs2sctw2*Qf*Qq*deng * ( 2.D0*  
     &    mq2*mf2*denzd + 2.D0*mq2*mf2*denz + p3p4*mq2*denzd + p3p4*mq2 
     &    *denz + p1p4*p2p3*denzd + p1p4*p2p3*denz + p1p3*p2p4*denzd +  
     &    p1p3*p2p4*denz + p1p2*mf2*denzd + p1p2*mf2*denz )             
      elmat2 = elmat2 + e2*ovall*gaf*gaq*gs2sctw2*Qf*Qq*deng * ( p1p4*  
     &    p2p3*denzd + p1p4*p2p3*denz - p1p3*p2p4*denzd - p1p3*p2p4*    
     &    denz )                                                        
                                                                        
                                                                        
 !***** number of lines = 43
