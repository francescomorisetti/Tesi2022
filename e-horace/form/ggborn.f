      elmat2 =                                                          
     &  + p2p3m2*e4 * ( 8.D0*p2p3*p2p4 )                                
      elmat2 = elmat2 + p1p3m2*e4 * ( 8.D0*p1p3*p1p4 )                  
      elmat2 = elmat2 + p2p3m1*p1p3m1*e4 * ( 16.D0*p2p3*p3p4 + 16.D0*   
     &    p1p3*p3p4 - 16.D0*p1p2*p3p4 )                                 
      elmat2 = elmat2 + mf4*p2p3m2*e4 * (  - 16.D0 )                    
      elmat2 = elmat2 + mf4*p1p3m2*e4 * (  - 16.D0 )                    
      elmat2 = elmat2 + mf4*p2p3m1*p1p3m1*e4 * (  - 8.D0 )              
      elmat2 = elmat2 + mf2*p2p3m2*e4 * (  - 8.D0*p3p4 + 8.D0*p2p4 + 16.
     &    D0*p2p3 )                                                     
      elmat2 = elmat2 + mf2*p1p3m2*e4 * (  - 8.D0*p3p4 + 8.D0*p1p4 + 16.
     &    D0*p1p3 )                                                     
      elmat2 = elmat2 + mf2*p2p3m1*p1p3m1*e4 * ( 8.D0*p3p4 - 8.D0*p2p4  
     &     + 16.D0*p2p3 - 8.D0*p1p4 + 16.D0*p1p3 - 8.D0*p1p2 )          
                                                                        
                                                                        
 !***** number of lines = 17
