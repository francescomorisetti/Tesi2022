      elmat2 =                                                          
     &  + Qw2*prop12m2*prop34m2 * ( 256.D0*p2p3*kp1*kp4 - 512.D0*p2p3*  
     &    p3p4*kp1 + 256.D0*p1p4*kp2*kp3 - 512.D0*p1p4*p3p4*kp2 + 256.D0
     &    *p1p4*p2p4*p3p4 + 256.D0*p1p4*p2p3*kp4 + 256.D0*p1p4*p2p3*kp3 
     &     - 256.D0*p1p4*p2p3*kp2 - 256.D0*p1p4*p2p3*kp1 - 256.D0*p1p4* 
     &    p2p3*p3p4 + 256.D0*p1p4*p2p3*p2p4 - 256.D0*p1p3*p2p4*kp4 -    
     &    256.D0*p1p3*p2p4*kp3 + 256.D0*p1p3*p2p4*kp2 + 256.D0*p1p3*    
     &    p2p4*kp1 + 256.D0*p1p3*p2p4*p3p4 + 256.D0*p1p3*p2p3*p3p4 +    
     &    256.D0*p1p3*p1p4*p2p3 + 256.D0*p1p2*p3p4*kp4 + 256.D0*p1p2*   
     &    p3p4*kp3 - 256.D0*p1p2*p3p4*kp2 - 256.D0*p1p2*p3p4*kp1 + 256.D
     &    0*p1p2*p2p4*p3p4 + 512.D0*p1p2*p2p3*kp4 - 512.D0*p1p2*p2p3*   
     &    p3p4 + 256.D0*p1p2*p2p3*p2p4 + 512.D0*p1p2*p1p4*kp3 - 512.D0* 
     &    p1p2*p1p4*p3p4 - 256.D0*p1p2*p1p4*p2p3 + 256.D0*p1p2*p1p3*    
     &    p3p4 + 256.D0*p1p2*p1p3*p2p4 + 256.D0*p1p2*p1p3*p1p4 - 256.D0 
     &    *p3p42*p1p2 - 256.D0*p2p42*p1p3 - 256.D0*p1p32*p2p4 - 256.D0* 
     &    p1p22*p3p4 )                                                  
      elmat2 = elmat2 + Qw2*md2*prop12m2*prop34m2 * (  - 128.D0*p3p4*   
     &    kp1 + 384.D0*p1p4*kp3 - 512.D0*p1p4*p3p4 - 128.D0*p1p4*p2p3   
     &     - 128.D0*p1p3*kp4 + 256.D0*p1p3*p3p4 + 128.D0*p1p3*p2p4 +    
     &    256.D0*p1p3*p1p4 - 256.D0*p1p2*p3p4 )                         
      elmat2 = elmat2 + Qw2*mup2*prop12m2*prop34m2 * (  - 128.D0*p3p4*  
     &    kp2 - 128.D0*p2p4*kp3 + 256.D0*p2p4*p3p4 + 384.D0*p2p3*kp4 -  
     &    512.D0*p2p3*p3p4 + 256.D0*p2p3*p2p4 - 128.D0*p1p4*p2p3 + 128.D
     &    0*p1p3*p2p4 - 256.D0*p1p2*p3p4 )                              
      elmat2 = elmat2 + Qw2*mup2*md2*prop12m2*prop34m2 * (  - 256.D0*   
     &    p3p4 )                                                        
      elmat2 = elmat2 + Qw2*me2*prop12m2*prop34m2 * (  - 384.D0*p2p3*   
     &    kp1 - 128.D0*p1p4*p2p3 + 128.D0*p1p3*kp2 + 128.D0*p1p3*p2p4   
     &     + 256.D0*p1p3*p2p3 + 128.D0*p1p2*kp3 - 256.D0*p1p2*p3p4 -    
     &    512.D0*p1p2*p2p3 + 256.D0*p1p2*p1p3 )                         
      elmat2 = elmat2 + Qw2*me2*md2*prop12m2*prop34m2 * ( 384.D0*p1p3 ) 
      elmat2 = elmat2 + Qw2*me2*mup2*prop12m2*prop34m2 * (  - 384.D0*   
     &    p2p3 )                                                        
      elmat2 = elmat2 + Qu2*kp1m2*mup2*prop34m2 * ( 256.D0*p2p3*kp4 -   
     &    256.D0*p1p4*p2p3 )                                            
      elmat2 = elmat2 + Qu2*kp1m1*prop34m2 * ( 256.D0*p2p3*kp4 )        
      elmat2 = elmat2 + Qd2*kp2m2*md2*prop34m2 * ( 256.D0*p1p4*kp3 -    
     &    256.D0*p1p4*p2p3 )                                            
      elmat2 = elmat2 + Qd2*kp2m1*prop34m2 * ( 256.D0*p1p4*kp3 )        
      elmat2 = elmat2 + Qe2*kp4m1*prop12m2 * ( 256.D0*p2p3*kp1 )        
      elmat2 = elmat2 + Qe2*me2*kp4m2*prop12m2 * (  - 256.D0*p2p3*kp1   
     &     - 256.D0*p1p4*p2p3 )                                         
      elmat2 = elmat2 + Qd*Qw*prop12d*prop34m2 * ( 256.D0*p1p4*kp3 -    
     &    128.D0*p1p4*p3p4 - 128.D0*p1p4*p2p3 + 128.D0*p1p3*p2p4 - 128.D
     &    0*p1p2*p3p4 )                                                 
      elmat2 = elmat2 + Qd*Qw*prop12*prop34m2 * ( 256.D0*p1p4*kp3 - 128.
     &    D0*p1p4*p3p4 - 128.D0*p1p4*p2p3 + 128.D0*p1p3*p2p4 - 128.D0*  
     &    p1p2*p3p4 )                                                   
      elmat2 = elmat2 + Qd*Qw*mup2*prop12d*prop34m2 * (  - 64.D0*p3p4 ) 
      elmat2 = elmat2 + Qd*Qw*mup2*prop12*prop34m2 * (  - 64.D0*p3p4 )  
      elmat2 = elmat2 + Qd*Qw*kp2m1*prop12d*prop34m2 * (  - 256.D0*p2p3 
     &    *p2p4*kp1 - 128.D0*p1p4*p2p4*kp3 + 128.D0*p1p4*p2p3*kp4 + 128.
     &    D0*p1p4*p2p3*kp3 - 128.D0*p1p3*p2p4*kp3 + 128.D0*p1p3*p2p3*   
     &    p2p4 + 128.D0*p1p2*p3p4*kp3 + 256.D0*p1p2*p2p3*kp4 - 128.D0*  
     &    p1p2*p2p3*p3p4 + 256.D0*p1p2*p1p4*kp3 - 256.D0*p1p2*p1p4*p2p3 
     &     - 128.D0*p2p32*p1p4 )                                        
      elmat2 = elmat2 + Qd*Qw*kp2m1*prop12*prop34m2 * (  - 256.D0*p2p3* 
     &    p2p4*kp1 - 128.D0*p1p4*p2p4*kp3 + 128.D0*p1p4*p2p3*kp4 + 128.D
     &    0*p1p4*p2p3*kp3 - 128.D0*p1p3*p2p4*kp3 + 128.D0*p1p3*p2p3*    
     &    p2p4 + 128.D0*p1p2*p3p4*kp3 + 256.D0*p1p2*p2p3*kp4 - 128.D0*  
     &    p1p2*p2p3*p3p4 + 256.D0*p1p2*p1p4*kp3 - 256.D0*p1p2*p1p4*p2p3 
     &     - 128.D0*p2p32*p1p4 )                                        
      elmat2 = elmat2 + Qd*Qw*kp2m1*md2*prop12d*prop34m2 * ( 64.D0*p3p4 
     &    *kp1 + 320.D0*p1p4*kp3 - 128.D0*p1p4*p3p4 - 192.D0*p1p4*p2p3  
     &     - 64.D0*p1p3*kp4 + 64.D0*p1p3*p2p4 + 128.D0*p1p3*p1p4 - 64.D0
     &    *p1p2*p3p4 )                                                  
      elmat2 = elmat2 + Qd*Qw*kp2m1*md2*prop12*prop34m2 * ( 64.D0*p3p4* 
     &    kp1 + 320.D0*p1p4*kp3 - 128.D0*p1p4*p3p4 - 192.D0*p1p4*p2p3   
     &     - 64.D0*p1p3*kp4 + 64.D0*p1p3*p2p4 + 128.D0*p1p3*p1p4 - 64.D0
     &    *p1p2*p3p4 )                                                  
      elmat2 = elmat2 + Qd*Qw*kp2m1*mup2*prop12d*prop34m2 * (  - 64.D0* 
     &    p2p4*kp3 + 64.D0*p2p3*kp4 )                                   
      elmat2 = elmat2 + Qd*Qw*kp2m1*mup2*prop12*prop34m2 * (  - 64.D0*  
     &    p2p4*kp3 + 64.D0*p2p3*kp4 )                                   
      elmat2 = elmat2 + Qd*Qw*kp2m1*mup2*md2*prop12d*prop34m2 * (  - 64.
     &    D0*p3p4 )                                                     
      elmat2 = elmat2 + Qd*Qw*kp2m1*mup2*md2*prop12*prop34m2 * (  - 64.D
     &    0*p3p4 )                                                      
      elmat2 = elmat2 + Qd*Qw*eps1234*prop12d*prop34m2 * ( 128.D0 )     
      elmat2 = elmat2 + Qd*Qw*eps1234*prop12*prop34m2 * (  - 128.D0 )   
      elmat2 = elmat2 + Qd*Qw*eps1234*kp2m1*prop12d*prop34m2 * (  - 128.
     &    D0*kp3 + 384.D0*p2p3 + 128.D0*p1p4 )                          
      elmat2 = elmat2 + Qd*Qw*eps1234*kp2m1*prop12*prop34m2 * ( 128.D0* 
     &    kp3 - 384.D0*p2p3 - 128.D0*p1p4 )                             
      elmat2 = elmat2 + Qd*Qw*eps1234*kp2m1*mup2*prop12d*prop34m2 * (   
     &    64.D0 )                                                       
      elmat2 = elmat2 + Qd*Qw*eps1234*kp2m1*mup2*prop12*prop34m2 * (    
     &     - 64.D0 )                                                    
      elmat2 = elmat2 + Qd*Qw*me2*kp2m1*prop12d*prop34m2 * ( 128.D0*    
     &    p1p2*kp3 - 128.D0*p1p2*p2p3 )                                 
      elmat2 = elmat2 + Qd*Qw*me2*kp2m1*prop12*prop34m2 * ( 128.D0*p1p2 
     &    *kp3 - 128.D0*p1p2*p2p3 )                                     
      elmat2 = elmat2 + Qd*Qw*me2*kp2m1*md2*prop12d*prop34m2 * ( 64.D0* 
     &    p1p3 )                                                        
      elmat2 = elmat2 + Qd*Qw*me2*kp2m1*md2*prop12*prop34m2 * ( 64.D0*  
     &    p1p3 )                                                        
      elmat2 = elmat2 + Qu*Qw*prop12d*prop34m2 * (  - 256.D0*p2p3*kp4   
     &     + 128.D0*p2p3*p3p4 + 128.D0*p1p4*p2p3 - 128.D0*p1p3*p2p4 +   
     &    128.D0*p1p2*p3p4 )                                            
      elmat2 = elmat2 + Qu*Qw*prop12*prop34m2 * (  - 256.D0*p2p3*kp4 +  
     &    128.D0*p2p3*p3p4 + 128.D0*p1p4*p2p3 - 128.D0*p1p3*p2p4 + 128.D
     &    0*p1p2*p3p4 )                                                 
      elmat2 = elmat2 + Qu*Qw*md2*prop12d*prop34m2 * ( 64.D0*p3p4 )     
      elmat2 = elmat2 + Qu*Qw*md2*prop12*prop34m2 * ( 64.D0*p3p4 )      
      elmat2 = elmat2 + Qu*Qw*kp1m1*prop12d*prop34m2 * (  - 128.D0*p1p4 
     &    *p2p3*kp4 - 128.D0*p1p4*p2p3*kp3 + 128.D0*p1p3*p2p4*kp4 + 128.
     &    D0*p1p3*p2p3*kp4 + 256.D0*p1p3*p1p4*kp2 - 128.D0*p1p3*p1p4*   
     &    p2p4 - 128.D0*p1p2*p3p4*kp4 - 256.D0*p1p2*p2p3*kp4 - 256.D0*  
     &    p1p2*p1p4*kp3 + 128.D0*p1p2*p1p4*p3p4 + 256.D0*p1p2*p1p4*p2p3 
     &     + 128.D0*p1p42*p2p3 )                                        
      elmat2 = elmat2 + Qu*Qw*kp1m1*prop12*prop34m2 * (  - 128.D0*p1p4* 
     &    p2p3*kp4 - 128.D0*p1p4*p2p3*kp3 + 128.D0*p1p3*p2p4*kp4 + 128.D
     &    0*p1p3*p2p3*kp4 + 256.D0*p1p3*p1p4*kp2 - 128.D0*p1p3*p1p4*    
     &    p2p4 - 128.D0*p1p2*p3p4*kp4 - 256.D0*p1p2*p2p3*kp4 - 256.D0*  
     &    p1p2*p1p4*kp3 + 128.D0*p1p2*p1p4*p3p4 + 256.D0*p1p2*p1p4*p2p3 
     &     + 128.D0*p1p42*p2p3 )                                        
      elmat2 = elmat2 + Qu*Qw*kp1m1*md2*prop12d*prop34m2 * (  - 64.D0*  
     &    p1p4*kp3 + 64.D0*p1p3*kp4 )                                   
      elmat2 = elmat2 + Qu*Qw*kp1m1*md2*prop12*prop34m2 * (  - 64.D0*   
     &    p1p4*kp3 + 64.D0*p1p3*kp4 )                                   
      elmat2 = elmat2 + Qu*Qw*kp1m1*mup2*prop12d*prop34m2 * (  - 64.D0* 
     &    p3p4*kp2 + 64.D0*p2p4*kp3 - 320.D0*p2p3*kp4 + 128.D0*p2p3*    
     &    p3p4 - 128.D0*p2p3*p2p4 + 192.D0*p1p4*p2p3 - 64.D0*p1p3*p2p4  
     &     + 64.D0*p1p2*p3p4 )                                          
      elmat2 = elmat2 + Qu*Qw*kp1m1*mup2*prop12*prop34m2 * (  - 64.D0*  
     &    p3p4*kp2 + 64.D0*p2p4*kp3 - 320.D0*p2p3*kp4 + 128.D0*p2p3*    
     &    p3p4 - 128.D0*p2p3*p2p4 + 192.D0*p1p4*p2p3 - 64.D0*p1p3*p2p4  
     &     + 64.D0*p1p2*p3p4 )                                          
      elmat2 = elmat2 + Qu*Qw*kp1m1*mup2*md2*prop12d*prop34m2 * ( 64.D0 
     &    *p3p4 )                                                       
      elmat2 = elmat2 + Qu*Qw*kp1m1*mup2*md2*prop12*prop34m2 * ( 64.D0* 
     &    p3p4 )                                                        
      elmat2 = elmat2 + Qu*Qw*eps1234*kp1m1*prop12d*prop34m2 * (  - 128.
     &    D0*kp4 + 128.D0*p2p3 + 256.D0*p1p4 - 128.D0*p1p3 + 128.D0*    
     &    p1p2 )                                                        
      elmat2 = elmat2 + Qu*Qw*eps1234*kp1m1*prop12*prop34m2 * ( 128.D0* 
     &    kp4 - 128.D0*p2p3 - 256.D0*p1p4 + 128.D0*p1p3 - 128.D0*p1p2 ) 
      elmat2 = elmat2 + Qu*Qw*eps1234*kp1m1*md2*prop12d*prop34m2 * ( 64.
     &    D0 )                                                          
      elmat2 = elmat2 + Qu*Qw*eps1234*kp1m1*md2*prop12*prop34m2 * (  -  
     &    64.D0 )                                                       
      elmat2 = elmat2 + Qu*Qw*eps1234*kp1m1*mup2*prop12d*prop34m2 * (   
     &    128.D0 )                                                      
      elmat2 = elmat2 + Qu*Qw*eps1234*kp1m1*mup2*prop12*prop34m2 * (    
     &     - 128.D0 )                                                   
      elmat2 = elmat2 + Qu*Qw*me2*prop12d*prop34m2 * ( 128.D0*p2p3 )    
      elmat2 = elmat2 + Qu*Qw*me2*prop12*prop34m2 * ( 128.D0*p2p3 )     
      elmat2 = elmat2 + Qu*Qw*me2*kp1m1*mup2*prop12d*prop34m2 * ( 64.D0 
     &    *p2p3 )                                                       
      elmat2 = elmat2 + Qu*Qw*me2*kp1m1*mup2*prop12*prop34m2 * ( 64.D0* 
     &    p2p3 )                                                        
      elmat2 = elmat2 + Qu*Qd*kp2m1*prop34m2 * ( 256.D0*p2p3*p2p4 - 256.
     &    D0*p1p4*p2p3 )                                                
      elmat2 = elmat2 + Qu*Qd*kp1m1*prop34m2 * (  - 256.D0*p1p4*p2p3 +  
     &    256.D0*p1p3*p1p4 )                                            
      elmat2 = elmat2 + Qu*Qd*kp1m1*kp2m1*prop34m2 * (  - 256.D0*p1p2*  
     &    p2p3*kp4 - 256.D0*p1p2*p1p4*kp3 + 512.D0*p1p2*p1p4*p2p3 )     
      elmat2 = elmat2 + Qe*Qw*prop34d*prop12m2 * (  - 256.D0*p2p3*kp1   
     &     - 128.D0*p1p4*p2p3 + 128.D0*p1p3*p2p4 - 128.D0*p1p2*p3p4 -   
     &    128.D0*p1p2*p2p3 )                                            
      elmat2 = elmat2 + Qe*Qw*prop34*prop12m2 * (  - 256.D0*p2p3*kp1 -  
     &    128.D0*p1p4*p2p3 + 128.D0*p1p3*p2p4 - 128.D0*p1p2*p3p4 - 128.D
     &    0*p1p2*p2p3 )                                                 
      elmat2 = elmat2 + Qe*Qw*mup2*prop34d*prop12m2 * (  - 128.D0*p2p3  
     &     )                                                            
      elmat2 = elmat2 + Qe*Qw*mup2*prop34*prop12m2 * (  - 128.D0*p2p3 ) 
      elmat2 = elmat2 + Qe*Qw*kp4m1*prop34d*prop12m2 * ( 256.D0*p2p3*   
     &    p3p4*kp1 - 128.D0*p2p3*p2p4*kp1 + 256.D0*p1p4*p3p4*kp2 - 256.D
     &    0*p1p4*p2p4*kp3 + 128.D0*p1p4*p2p3*kp2 + 128.D0*p1p4*p2p3*kp1 
     &     + 256.D0*p1p4*p2p3*p3p4 - 128.D0*p1p3*p2p4*kp1 - 128.D0*p1p3 
     &    *p1p4*p2p4 + 128.D0*p1p2*p3p4*kp1 + 128.D0*p1p2*p1p4*p3p4 +   
     &    128.D0*p1p42*p2p3 )                                           
      elmat2 = elmat2 + Qe*Qw*kp4m1*prop34*prop12m2 * ( 256.D0*p2p3*    
     &    p3p4*kp1 - 128.D0*p2p3*p2p4*kp1 + 256.D0*p1p4*p3p4*kp2 - 256.D
     &    0*p1p4*p2p4*kp3 + 128.D0*p1p4*p2p3*kp2 + 128.D0*p1p4*p2p3*kp1 
     &     + 256.D0*p1p4*p2p3*p3p4 - 128.D0*p1p3*p2p4*kp1 - 128.D0*p1p3 
     &    *p1p4*p2p4 + 128.D0*p1p2*p3p4*kp1 + 128.D0*p1p2*p1p4*p3p4 +   
     &    128.D0*p1p42*p2p3 )                                           
      elmat2 = elmat2 + Qe*Qw*kp4m1*md2*prop34d*prop12m2 * ( 128.D0*    
     &    p3p4*kp1 + 128.D0*p1p4*p3p4 )                                 
      elmat2 = elmat2 + Qe*Qw*kp4m1*md2*prop34*prop12m2 * ( 128.D0*p3p4 
     &    *kp1 + 128.D0*p1p4*p3p4 )                                     
      elmat2 = elmat2 + Qe*Qw*eps1234*kp4m1*prop34d*prop12m2 * ( 128.D0 
     &    *kp1 + 128.D0*p3p4 - 128.D0*p2p4 + 128.D0*p2p3 + 256.D0*p1p4  
     &     )                                                            
      elmat2 = elmat2 + Qe*Qw*eps1234*kp4m1*prop34*prop12m2 * (  - 128.D
     &    0*kp1 - 128.D0*p3p4 + 128.D0*p2p4 - 128.D0*p2p3 - 256.D0*p1p4 
     &     )                                                            
      elmat2 = elmat2 + Qe*Qw*me2*kp4m1*prop34d*prop12m2 * ( 320.D0*    
     &    p2p3*kp1 + 192.D0*p1p4*p2p3 - 64.D0*p1p3*kp2 - 64.D0*p1p3*    
     &    p2p4 - 128.D0*p1p3*p2p3 + 64.D0*p1p2*kp3 + 64.D0*p1p2*p3p4 +  
     &    128.D0*p1p2*p2p3 )                                            
      elmat2 = elmat2 + Qe*Qw*me2*kp4m1*prop34*prop12m2 * ( 320.D0*p2p3 
     &    *kp1 + 192.D0*p1p4*p2p3 - 64.D0*p1p3*kp2 - 64.D0*p1p3*p2p4 -  
     &    128.D0*p1p3*p2p3 + 64.D0*p1p2*kp3 + 64.D0*p1p2*p3p4 + 128.D0* 
     &    p1p2*p2p3 )                                                   
      elmat2 = elmat2 + Qe*Qw*me2*kp4m1*md2*prop34d*prop12m2 * (  - 64.D
     &    0*p1p3 )                                                      
      elmat2 = elmat2 + Qe*Qw*me2*kp4m1*md2*prop34*prop12m2 * (  - 64.D0
     &    *p1p3 )                                                       
      elmat2 = elmat2 + Qe*Qw*me2*kp4m1*mup2*prop34d*prop12m2 * ( 64.D0 
     &    *p2p3 )                                                       
      elmat2 = elmat2 + Qe*Qw*me2*kp4m1*mup2*prop34*prop12m2 * ( 64.D0* 
     &    p2p3 )                                                        
      elmat2 = elmat2 + Qe*Qw*me2*eps1234*kp4m1*prop34d*prop12m2 * (    
     &    128.D0 )                                                      
      elmat2 = elmat2 + Qe*Qw*me2*eps1234*kp4m1*prop34*prop12m2 * (  -  
     &    128.D0 )                                                      
      elmat2 = elmat2 + Qe*Qd*kp4m1*prop12d*prop34 * ( 128.D0*p1p4*p3p4 
     &     + 128.D0*p1p4*p2p3 )                                         
      elmat2 = elmat2 + Qe*Qd*kp4m1*prop12*prop34d * ( 128.D0*p1p4*p3p4 
     &     + 128.D0*p1p4*p2p3 )                                         
      elmat2 = elmat2 + Qe*Qd*kp2m1*prop12d*prop34 * (  - 128.D0*p1p4*  
     &    p2p3 - 128.D0*p1p2*p2p3 )                                     
      elmat2 = elmat2 + Qe*Qd*kp2m1*prop12*prop34d * (  - 128.D0*p1p4*  
     &    p2p3 - 128.D0*p1p2*p2p3 )                                     
      elmat2 = elmat2 + Qe*Qd*kp2m1*kp4m1*prop12d*prop34 * ( 128.D0*    
     &    p2p3*p2p4*kp1 - 128.D0*p1p4*p2p4*kp3 + 256.D0*p1p4*p2p3*p2p4  
     &     )                                                            
      elmat2 = elmat2 + Qe*Qd*kp2m1*kp4m1*prop12*prop34d * ( 128.D0*    
     &    p2p3*p2p4*kp1 - 128.D0*p1p4*p2p4*kp3 + 256.D0*p1p4*p2p3*p2p4  
     &     )                                                            
      elmat2 = elmat2 + Qe*Qd*eps1234*kp2m1*prop12d*prop34 * (  - 128.D0
     &     )                                                            
      elmat2 = elmat2 + Qe*Qd*eps1234*kp2m1*kp4m1*prop12d*prop34 * (    
     &     - 128.D0*p3p4 + 128.D0*p2p4 - 128.D0*p2p3 )                  
      elmat2 = elmat2 + Qe*Qd*eps1234*kp2m1*kp4m1*prop12*prop34d * (    
     &    128.D0*p2p3 + 128.D0*p1p4 )                                   
      elmat2 = elmat2 + Qe*Qd*me2*eps1234*kp2m1*kp4m1*prop12d*prop34    
     &  * (  - 128.D0 )                                                 
      elmat2 = elmat2 + Qe*Qu*prop12d*prop34 * ( 256.D0*p2p3 )          
      elmat2 = elmat2 + Qe*Qu*prop12*prop34d * ( 256.D0*p2p3 )          
      elmat2 = elmat2 + Qe*Qu*kp4m1*prop12d*prop34 * (  - 256.D0*p1p4*  
     &    p2p3 )                                                        
      elmat2 = elmat2 + Qe*Qu*kp4m1*prop12*prop34d * (  - 256.D0*p1p4*  
     &    p2p3 )                                                        
      elmat2 = elmat2 + Qe*Qu*kp1m1*prop12d*prop34 * ( 256.D0*p1p4*p2p3 
     &     )                                                            
      elmat2 = elmat2 + Qe*Qu*kp1m1*prop12*prop34d * ( 256.D0*p1p4*p2p3 
     &     )                                                            
      elmat2 = elmat2 + Qe*Qu*kp1m1*mup2*prop12d*prop34 * ( 128.D0*p2p3 
     &     )                                                            
      elmat2 = elmat2 + Qe*Qu*kp1m1*mup2*prop12*prop34d * ( 128.D0*p2p3 
     &     )                                                            
      elmat2 = elmat2 + Qe*Qu*kp1m1*kp4m1*prop12d*prop34 * (  - 256.D0* 
     &    p1p42*p2p3 )                                                  
      elmat2 = elmat2 + Qe*Qu*kp1m1*kp4m1*prop12*prop34d * (  - 256.D0* 
     &    p1p42*p2p3 )                                                  
      elmat2 = elmat2 + Qe*Qu*me2*kp4m1*prop12d*prop34 * (  - 128.D0*   
     &    p2p3 )                                                        
      elmat2 = elmat2 + Qe*Qu*me2*kp4m1*prop12*prop34d * (  - 128.D0*   
     &    p2p3 )                                                        
                                                                        
                                                                        
 !***** number of lines = 255
