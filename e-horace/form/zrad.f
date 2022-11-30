      elmat2 =                                                          
     &  + gvf*gvq*gs2sctw2*Qf*Qq*mq4*mf2*Qq2*e4*dg34*dz34m2*s34mmz2     
     &  * (  - 128.D0*kp2m2 - 128.D0*kp1m2 )                            
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*Qf*Qq*mq2*mf4*Qf2*e4*dg12*     
     & dz12m2*s12mmz2 * (  - 128.D0*kp4m2 - 128.D0*kp3m2 )              
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*Qf*Qq*mq2*mf2*Qq2*e4*dg34*     
     & dz34m2*s34mmz2 * ( 64.D0*kp2*kp1m2 + 64.D0*kp1*kp2m2 )           
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*Qf*Qq*mq2*mf2*Qf2*e4*dg12*     
     & dz12m2*s12mmz2 * (  - 64.D0*kp4*kp3m2 - 64.D0*kp3*kp4m2 )        
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg34*     
     & dz12m2*s12mmz2 * ( 32.D0*kp2 - 32.D0*kp1 )                       
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg12*     
     & dz34m2*s34mmz2 * ( 32.D0*kp2 - 32.D0*kp1 )                       
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp4m1*Qf*Qq*mq2*Qf2*e4*dg12*   
     & dz12m2*s12mmz2 * ( 64.D0*kp3 )                                   
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp4m1*Qf*Qq*mq2*mf2*Qf2*e4*    
     & dg12*dz12m2*s12mmz2 * (  - 64.D0 )                               
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg34*     
     & dz12m2*s12mmz2 * (  - 32.D0*kp2 + 32.D0*kp1 )                    
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg12*     
     & dz34m2*s34mmz2 * (  - 32.D0*kp2 + 32.D0*kp1 )                    
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp3m1*Qf*Qq*mq2*Qf2*e4*dg12*   
     & dz12m2*s12mmz2 * ( 64.D0*kp4 )                                   
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp3m1*Qf*Qq*mq2*mf2*Qf2*e4*    
     & dg12*dz12m2*s12mmz2 * (  - 64.D0 )                               
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*mf2*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * (  - 128.D0*kp1*kp2 )                      
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*mq2*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * ( 128.D0*p3p42 )                           
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg34*     
     & dz12m2*s12mmz2 * ( 32.D0*kp4 - 32.D0*kp3 )                       
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg12*     
     & dz34m2*s34mmz2 * ( 32.D0*kp4 - 32.D0*kp3 )                       
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp2m1*Qf*Qq*mf2*Qq2*e4*dg34*   
     & dz34m2*s34mmz2 * ( 64.D0*kp1 )                                   
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp2m1*Qf*Qq*mq2*mf2*Qq2*e4*    
     & dg34*dz34m2*s34mmz2 * ( 64.D0 )                                  
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg34*     
     & dz12m2*s12mmz2 * (  - 32.D0*kp4 + 32.D0*kp3 )                    
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg12*     
     & dz34m2*s34mmz2 * (  - 32.D0*kp4 + 32.D0*kp3 )                    
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp1m1*Qf*Qq*mf2*Qq2*e4*dg34*   
     & dz34m2*s34mmz2 * ( 64.D0*kp2 )                                   
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp1m1*Qf*Qq*mq2*mf2*Qq2*e4*    
     & dg34*dz34m2*s34mmz2 * ( 64.D0 )                                  
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*mf2*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * ( 128.D0*p1p22 )                           
      elmat2 = elmat2 + gvf*gvq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*mq2*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * (  - 128.D0*kp3*kp4 )                      
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg34*     
     & dz12m2*s12mmz2 * ( 32.D0*kp2 + 32.D0*kp1 )                       
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg12*     
     & dz34m2*s34mmz2 * (  - 32.D0*kp2 - 32.D0*kp1 )                    
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg34*     
     & dz12m2*s12mmz2 * ( 32.D0*kp2 + 32.D0*kp1 )                       
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg12*     
     & dz34m2*s34mmz2 * (  - 32.D0*kp2 - 32.D0*kp1 )                    
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg34*     
     & dz12m2*s12mmz2 * (  - 32.D0*kp4 - 32.D0*kp3 )                    
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg12*     
     & dz34m2*s34mmz2 * ( 32.D0*kp4 + 32.D0*kp3 )                       
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp2m1*kp4m1*mq2*mf2*Qf2*Qq2*e4 
     & *dg34*dz12m2*s12mmz2 * ( 32.D0*kp3 + 32.D0*kp1 )                 
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp2m1*kp4m1*mq2*mf2*Qf2*Qq2*e4 
     & *dg12*dz34m2*s34mmz2 * (  - 32.D0*kp3 - 32.D0*kp1 )              
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp2m1*kp3m1*mq2*mf2*Qf2*Qq2*e4 
     & *dg34*dz12m2*s12mmz2 * ( 32.D0*kp4 + 32.D0*kp1 )                 
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp2m1*kp3m1*mq2*mf2*Qf2*Qq2*e4 
     & *dg12*dz34m2*s34mmz2 * (  - 32.D0*kp4 - 32.D0*kp1 )              
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg34*     
     & dz12m2*s12mmz2 * (  - 32.D0*kp4 - 32.D0*kp3 )                    
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg12*     
     & dz34m2*s34mmz2 * ( 32.D0*kp4 + 32.D0*kp3 )                       
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp1m1*kp4m1*mq2*mf2*Qf2*Qq2*e4 
     & *dg34*dz12m2*s12mmz2 * ( 32.D0*kp3 + 32.D0*kp2 )                 
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp1m1*kp4m1*mq2*mf2*Qf2*Qq2*e4 
     & *dg12*dz34m2*s34mmz2 * (  - 32.D0*kp3 - 32.D0*kp2 )              
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp1m1*kp3m1*mq2*mf2*Qf2*Qq2*e4 
     & *dg34*dz12m2*s12mmz2 * ( 32.D0*kp4 + 32.D0*kp2 )                 
      elmat2 = elmat2 + gaf*gaq*gs2sctw2*kp1m1*kp3m1*mq2*mf2*Qf2*Qq2*e4 
     & *dg12*dz34m2*s34mmz2 * (  - 32.D0*kp4 - 32.D0*kp2 )              
      elmat2 = elmat2 + p3p4*gvf*gvq*gs2sctw2*Qf*Qq*mq4*Qq2*e4*dg34*    
     & dz34m2*s34mmz2 * (  - 64.D0*kp2m2 - 64.D0*kp1m2 )                
      elmat2 = elmat2 + p3p4*gvf*gvq*gs2sctw2*Qf*Qq*mq2*mf2*Qf2*e4*dg12 
     & *dz12m2*s12mmz2 * (  - 64.D0*kp4m2 - 64.D0*kp3m2 )               
      elmat2 = elmat2 + p3p4*gvf*gvq*gs2sctw2*kp4m1*Qf*Qq*mq2*Qf2*e4*   
     & dg12*dz12m2*s12mmz2 * ( 128.D0 )                                 
      elmat2 = elmat2 + p3p4*gvf*gvq*gs2sctw2*kp3m1*Qf*Qq*mq2*Qf2*e4*   
     & dg12*dz12m2*s12mmz2 * ( 128.D0 )                                 
      elmat2 = elmat2 + p3p4*gvf*gvq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*mq2*mf2 
     & *Qf2*e4*dg12*dz12m2*s12mmz2 * ( 256.D0 )                         
      elmat2 = elmat2 + p3p4*gaf*gaq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + p3p4*gaf*gaq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + p3p4*gaf*gaq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + p3p4*gaf*gaq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*Qf2*Qq2*e4*dg34*dz12m2*   
     & s12mmz2 * (  - 64.D0 )                                           
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*Qf2*Qq2*e4*dg12*dz34m2*   
     & s34mmz2 * (  - 64.D0 )                                           
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*dg12*    
     & dz12m2*s12mmz2 * (  - 64.D0*kp1*kp3m2 )                          
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*dg34*    
     & dz34m2*s34mmz2 * ( 64.D0*kp3*kp1m2 )                             
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp3m1*Qf*Qq*Qf2*e4*dg12*  
     & dz12m2*s12mmz2 * ( 64.D0*kp1 )                                   
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*mf2*Qf2*Qq2*  
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp1 )                        
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*mf2*Qf2*Qq2*  
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp1 )                        
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*mq2*Qf2*Qq2*  
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp3 )                           
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*mq2*Qf2*Qq2*  
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp3 )                           
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*mq2*mf2*Qf2*  
     & Qq2*e4*dg34*dz12m2*s12mmz2 * ( 128.D0 )                          
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*mq2*mf2*Qf2*  
     & Qq2*e4*dg12*dz34m2*s34mmz2 * ( 128.D0 )                          
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp1m1*Qf*Qq*Qq2*e4*dg34*  
     & dz34m2*s34mmz2 * ( 64.D0*kp3 )                                   
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2*e4*   
     & dg34*dz12m2*s12mmz2 * ( 64.D0*p1p32 )                            
      elmat2 = elmat2 + p2p4*gvf*gvq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2*e4*   
     & dg12*dz34m2*s34mmz2 * ( 64.D0*p1p32 )                            
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*Qf2*Qq2*e4*dg34*dz12m2*   
     & s12mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*Qf2*Qq2*e4*dg12*dz34m2*   
     & s34mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*dg12*    
     & dz12m2*s12mmz2 * ( 64.D0*kp1*kp3m2 )                             
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*dg34*    
     & dz34m2*s34mmz2 * (  - 64.D0*kp3*kp1m2 )                          
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*kp3m1*Qf*Qq*Qf2*e4*dg12*  
     & dz12m2*s12mmz2 * (  - 64.D0*kp1 )                                
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*kp1m1*Qf*Qq*Qq2*e4*dg34*  
     & dz34m2*s34mmz2 * (  - 64.D0*kp3 )                                
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2*e4*   
     & dg34*dz12m2*s12mmz2 * (  - 64.D0*p1p32 )                         
      elmat2 = elmat2 + p2p4*gaf*gaq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2*e4*   
     & dg12*dz34m2*s34mmz2 * (  - 64.D0*p1p32 )                         
      elmat2 = elmat2 + p2p4*p3p4*gvf*gvq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p2p4*p3p4*gvf*gvq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p2p4*p3p4*gvf*gvq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*   
     & Qf2*e4*dg12*dz12m2*s12mmz2 * ( 64.D0*kp1 )                       
      elmat2 = elmat2 + p2p4*p3p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*mq2*Qf2* 
     & Qq2*e4*dg34*dz12m2*s12mmz2 * ( 64.D0 )                           
      elmat2 = elmat2 + p2p4*p3p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*mq2*Qf2* 
     & Qq2*e4*dg12*dz34m2*s34mmz2 * ( 64.D0 )                           
      elmat2 = elmat2 + p2p4*p3p4*gaf*gaq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p2p4*p3p4*gaf*gaq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p2p4*p3p4*gaf*gaq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*   
     & Qf2*e4*dg12*dz12m2*s12mmz2 * (  - 64.D0*kp1 )                    
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*Qf2*Qq2*e4*dg34*dz12m2*   
     & s12mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*Qf2*Qq2*e4*dg12*dz34m2*   
     & s34mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*dg12*    
     & dz12m2*s12mmz2 * (  - 64.D0*kp1*kp4m2 )                          
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*dg34*    
     & dz34m2*s34mmz2 * ( 64.D0*kp4*kp1m2 )                             
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp4m1*Qf*Qq*Qf2*e4*dg12*  
     & dz12m2*s12mmz2 * ( 64.D0*kp1 )                                   
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp2m1*kp3m1*mf2*Qf2*Qq2*  
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp1 )                           
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp2m1*kp3m1*mf2*Qf2*Qq2*  
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp1 )                           
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp2m1*kp3m1*mq2*Qf2*Qq2*  
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp4 )                        
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp2m1*kp3m1*mq2*Qf2*Qq2*  
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp4 )                        
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp2m1*kp3m1*mq2*mf2*Qf2*  
     & Qq2*e4*dg34*dz12m2*s12mmz2 * (  - 128.D0 )                       
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp2m1*kp3m1*mq2*mf2*Qf2*  
     & Qq2*e4*dg12*dz34m2*s34mmz2 * (  - 128.D0 )                       
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp1m1*Qf*Qq*Qq2*e4*dg34*  
     & dz34m2*s34mmz2 * ( 64.D0*kp4 )                                   
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2*e4*   
     & dg34*dz12m2*s12mmz2 * (  - 64.D0*p1p42 )                         
      elmat2 = elmat2 + p2p3*gvf*gvq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2*e4*   
     & dg12*dz34m2*s34mmz2 * (  - 64.D0*p1p42 )                         
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*Qf2*Qq2*e4*dg34*dz12m2*   
     & s12mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*Qf2*Qq2*e4*dg12*dz34m2*   
     & s34mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*dg12*    
     & dz12m2*s12mmz2 * (  - 64.D0*kp1*kp4m2 )                          
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*dg34*    
     & dz34m2*s34mmz2 * ( 64.D0*kp4*kp1m2 )                             
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*kp4m1*Qf*Qq*Qf2*e4*dg12*  
     & dz12m2*s12mmz2 * ( 64.D0*kp1 )                                   
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*kp1m1*Qf*Qq*Qq2*e4*dg34*  
     & dz34m2*s34mmz2 * ( 64.D0*kp4 )                                   
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2*e4*   
     & dg34*dz12m2*s12mmz2 * (  - 64.D0*p1p42 )                         
      elmat2 = elmat2 + p2p3*gaf*gaq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2*e4*   
     & dg12*dz34m2*s34mmz2 * (  - 64.D0*p1p42 )                         
      elmat2 = elmat2 + p2p3*p3p4*gvf*gvq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p2p3*p3p4*gvf*gvq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p2p3*p3p4*gvf*gvq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*   
     & Qf2*e4*dg12*dz12m2*s12mmz2 * ( 64.D0*kp1 )                       
      elmat2 = elmat2 + p2p3*p3p4*gvf*gvq*gs2sctw2*kp2m1*kp3m1*mq2*Qf2* 
     & Qq2*e4*dg34*dz12m2*s12mmz2 * (  - 64.D0 )                        
      elmat2 = elmat2 + p2p3*p3p4*gvf*gvq*gs2sctw2*kp2m1*kp3m1*mq2*Qf2* 
     & Qq2*e4*dg12*dz34m2*s34mmz2 * (  - 64.D0 )                        
      elmat2 = elmat2 + p2p3*p3p4*gaf*gaq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p2p3*p3p4*gaf*gaq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p2p3*p3p4*gaf*gaq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*   
     & Qf2*e4*dg12*dz12m2*s12mmz2 * ( 64.D0*kp1 )                       
      elmat2 = elmat2 + p2p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*Qf*Qq*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * ( 128.D0 )                                 
      elmat2 = elmat2 + p2p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp1 )                           
      elmat2 = elmat2 + p2p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp1 )                           
      elmat2 = elmat2 + p2p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp1 )                        
      elmat2 = elmat2 + p2p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp1 )                        
      elmat2 = elmat2 + p2p3*p2p4*gaf*gaq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp1 )                           
      elmat2 = elmat2 + p2p3*p2p4*gaf*gaq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp1 )                           
      elmat2 = elmat2 + p2p3*p2p4*gaf*gaq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp1 )                           
      elmat2 = elmat2 + p2p3*p2p4*gaf*gaq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp1 )                           
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*Qf2*Qq2*e4*dg34*dz12m2*   
     & s12mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*Qf2*Qq2*e4*dg12*dz34m2*   
     & s34mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*dg12*    
     & dz12m2*s12mmz2 * (  - 64.D0*kp2*kp3m2 )                          
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*dg34*    
     & dz34m2*s34mmz2 * ( 64.D0*kp3*kp2m2 )                             
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp3m1*Qf*Qq*Qf2*e4*dg12*  
     & dz12m2*s12mmz2 * ( 64.D0*kp2 )                                   
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp2m1*Qf*Qq*Qq2*e4*dg34*  
     & dz34m2*s34mmz2 * ( 64.D0*kp3 )                                   
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2*e4*   
     & dg34*dz12m2*s12mmz2 * (  - 64.D0*p2p32 )                         
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2*e4*   
     & dg12*dz34m2*s34mmz2 * (  - 64.D0*p2p32 )                         
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*mf2*Qf2*Qq2*  
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp2 )                           
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*mf2*Qf2*Qq2*  
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp2 )                           
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*mq2*Qf2*Qq2*  
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp3 )                        
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*mq2*Qf2*Qq2*  
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp3 )                        
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*mq2*mf2*Qf2*  
     & Qq2*e4*dg34*dz12m2*s12mmz2 * (  - 128.D0 )                       
      elmat2 = elmat2 + p1p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*mq2*mf2*Qf2*  
     & Qq2*e4*dg12*dz34m2*s34mmz2 * (  - 128.D0 )                       
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*Qf2*Qq2*e4*dg34*dz12m2*   
     & s12mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*Qf2*Qq2*e4*dg12*dz34m2*   
     & s34mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*dg12*    
     & dz12m2*s12mmz2 * (  - 64.D0*kp2*kp3m2 )                          
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*dg34*    
     & dz34m2*s34mmz2 * ( 64.D0*kp3*kp2m2 )                             
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*kp3m1*Qf*Qq*Qf2*e4*dg12*  
     & dz12m2*s12mmz2 * ( 64.D0*kp2 )                                   
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*kp2m1*Qf*Qq*Qq2*e4*dg34*  
     & dz34m2*s34mmz2 * ( 64.D0*kp3 )                                   
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2*e4*   
     & dg34*dz12m2*s12mmz2 * (  - 64.D0*p2p32 )                         
      elmat2 = elmat2 + p1p4*gaf*gaq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2*e4*   
     & dg12*dz34m2*s34mmz2 * (  - 64.D0*p2p32 )                         
      elmat2 = elmat2 + p1p4*p3p4*gvf*gvq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p3p4*gvf*gvq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p3p4*gvf*gvq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*   
     & Qf2*e4*dg12*dz12m2*s12mmz2 * ( 64.D0*kp2 )                       
      elmat2 = elmat2 + p1p4*p3p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*mq2*Qf2* 
     & Qq2*e4*dg34*dz12m2*s12mmz2 * (  - 64.D0 )                        
      elmat2 = elmat2 + p1p4*p3p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*mq2*Qf2* 
     & Qq2*e4*dg12*dz34m2*s34mmz2 * (  - 64.D0 )                        
      elmat2 = elmat2 + p1p4*p3p4*gaf*gaq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p3p4*gaf*gaq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p3p4*gaf*gaq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*   
     & Qf2*e4*dg12*dz12m2*s12mmz2 * ( 64.D0*kp2 )                       
      elmat2 = elmat2 + p1p4*p2p4*gvf*gvq*gs2sctw2*kp4m1*Qf*Qq*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * (  - 128.D0 )                              
      elmat2 = elmat2 + p1p4*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp3 )                        
      elmat2 = elmat2 + p1p4*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp3 )                        
      elmat2 = elmat2 + p1p4*p2p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp3 )                           
      elmat2 = elmat2 + p1p4*p2p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp3 )                           
      elmat2 = elmat2 + p1p4*p2p4*gaf*gaq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp3 )                        
      elmat2 = elmat2 + p1p4*p2p4*gaf*gaq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp3 )                        
      elmat2 = elmat2 + p1p4*p2p4*gaf*gaq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp3 )                        
      elmat2 = elmat2 + p1p4*p2p4*gaf*gaq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp3 )                        
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*    
     & dg12*dz12m2*s12mmz2 * (  - 64.D0*kp4m2 - 64.D0*kp3m2 )           
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*    
     & dg34*dz34m2*s34mmz2 * (  - 64.D0*kp2m2 - 64.D0*kp1m2 )           
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp4m1*Qf*Qq*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * ( 64.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp3m1*Qf*Qq*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * ( 64.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp2m1*Qf*Qq*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * (  - 64.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gvf*gvq*gs2sctw2*kp1m1*Qf*Qq*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * (  - 64.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*    
     & dg12*dz12m2*s12mmz2 * (  - 64.D0*kp4m2 - 64.D0*kp3m2 )           
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*    
     & dg34*dz34m2*s34mmz2 * (  - 64.D0*kp2m2 - 64.D0*kp1m2 )           
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp4m1*Qf*Qq*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * ( 64.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp3m1*Qf*Qq*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * ( 64.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp2m1*Qf*Qq*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * (  - 64.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p4*p2p3*gaf*gaq*gs2sctw2*kp1m1*Qf*Qq*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * (  - 64.D0 )                               
      elmat2 = elmat2 + p1p4*p2p3*p3p4*gvf*gvq*gs2sctw2*kp3m1*kp4m1*Qf* 
     & Qq*Qf2*e4*dg12*dz12m2*s12mmz2 * ( 128.D0 )                       
      elmat2 = elmat2 + p1p4*p2p3*p3p4*gaf*gaq*gs2sctw2*kp3m1*kp4m1*Qf* 
     & Qq*Qf2*e4*dg12*dz12m2*s12mmz2 * ( 128.D0 )                       
      elmat2 = elmat2 + p1p4*p2p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*Qf2 
     & *Qq2*e4*dg34*dz12m2*s12mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p4*p2p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*Qf2 
     & *Qq2*e4*dg12*dz34m2*s34mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p4*p2p3*p2p4*gaf*gaq*gs2sctw2*kp2m1*kp4m1*Qf2 
     & *Qq2*e4*dg34*dz12m2*s12mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p4*p2p3*p2p4*gaf*gaq*gs2sctw2*kp2m1*kp4m1*Qf2 
     & *Qq2*e4*dg12*dz34m2*s34mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*Qf2*Qq2*e4*dg34*dz12m2*   
     & s12mmz2 * (  - 64.D0 )                                           
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*Qf2*Qq2*e4*dg12*dz34m2*   
     & s34mmz2 * (  - 64.D0 )                                           
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*dg12*    
     & dz12m2*s12mmz2 * (  - 64.D0*kp2*kp4m2 )                          
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*dg34*    
     & dz34m2*s34mmz2 * ( 64.D0*kp4*kp2m2 )                             
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp4m1*Qf*Qq*Qf2*e4*dg12*  
     & dz12m2*s12mmz2 * ( 64.D0*kp2 )                                   
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp2m1*Qf*Qq*Qq2*e4*dg34*  
     & dz34m2*s34mmz2 * ( 64.D0*kp4 )                                   
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2*e4*   
     & dg34*dz12m2*s12mmz2 * ( 64.D0*p2p42 )                            
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2*e4*   
     & dg12*dz34m2*s34mmz2 * ( 64.D0*p2p42 )                            
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp1m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*mf2*Qf2*Qq2*  
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp2 )                        
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*mf2*Qf2*Qq2*  
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp2 )                        
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*mq2*Qf2*Qq2*  
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp4 )                           
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*mq2*Qf2*Qq2*  
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp4 )                           
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*mq2*mf2*Qf2*  
     & Qq2*e4*dg34*dz12m2*s12mmz2 * ( 128.D0 )                          
      elmat2 = elmat2 + p1p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*mq2*mf2*Qf2*  
     & Qq2*e4*dg12*dz34m2*s34mmz2 * ( 128.D0 )                          
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*Qf2*Qq2*e4*dg34*dz12m2*   
     & s12mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*Qf2*Qq2*e4*dg12*dz34m2*   
     & s34mmz2 * ( 64.D0 )                                              
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*dg12*    
     & dz12m2*s12mmz2 * ( 64.D0*kp2*kp4m2 )                             
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*dg34*    
     & dz34m2*s34mmz2 * (  - 64.D0*kp4*kp2m2 )                          
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * (  - 32.D0 )                                   
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*kp4m1*Qf*Qq*Qf2*e4*dg12*  
     & dz12m2*s12mmz2 * (  - 64.D0*kp2 )                                
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*kp2m1*mq2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 32.D0 )                                      
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*kp2m1*Qf*Qq*Qq2*e4*dg34*  
     & dz34m2*s34mmz2 * (  - 64.D0*kp4 )                                
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2*e4*   
     & dg34*dz12m2*s12mmz2 * (  - 64.D0*p2p42 )                         
      elmat2 = elmat2 + p1p3*gaf*gaq*gs2sctw2*kp2m1*kp4m1*Qf2*Qq2*e4*   
     & dg12*dz34m2*s34mmz2 * (  - 64.D0*p2p42 )                         
      elmat2 = elmat2 + p1p3*p3p4*gvf*gvq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p3*p3p4*gvf*gvq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p3*p3p4*gvf*gvq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*   
     & Qf2*e4*dg12*dz12m2*s12mmz2 * ( 64.D0*kp2 )                       
      elmat2 = elmat2 + p1p3*p3p4*gvf*gvq*gs2sctw2*kp1m1*kp3m1*mq2*Qf2* 
     & Qq2*e4*dg34*dz12m2*s12mmz2 * ( 64.D0 )                           
      elmat2 = elmat2 + p1p3*p3p4*gvf*gvq*gs2sctw2*kp1m1*kp3m1*mq2*Qf2* 
     & Qq2*e4*dg12*dz34m2*s34mmz2 * ( 64.D0 )                           
      elmat2 = elmat2 + p1p3*p3p4*gaf*gaq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p3*p3p4*gaf*gaq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p3*p3p4*gaf*gaq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*   
     & Qf2*e4*dg12*dz12m2*s12mmz2 * (  - 64.D0*kp2 )                    
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*    
     & dg12*dz12m2*s12mmz2 * (  - 64.D0*kp4m2 - 64.D0*kp3m2 )           
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*    
     & dg34*dz34m2*s34mmz2 * (  - 64.D0*kp2m2 - 64.D0*kp1m2 )           
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp4m1*Qf*Qq*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * ( 64.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp3m1*Qf*Qq*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * ( 64.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*Qf*Qq*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * (  - 64.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gvf*gvq*gs2sctw2*kp1m1*Qf*Qq*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * (  - 64.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*Qf*Qq*mf2*Qf2*e4*    
     & dg12*dz12m2*s12mmz2 * ( 64.D0*kp4m2 + 64.D0*kp3m2 )              
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*Qf*Qq*mq2*Qq2*e4*    
     & dg34*dz34m2*s34mmz2 * ( 64.D0*kp2m2 + 64.D0*kp1m2 )              
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp4m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp4m1*Qf*Qq*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * (  - 64.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp3m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp3m1*Qf*Qq*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * (  - 64.D0 )                               
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp2m1*Qf*Qq*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * ( 64.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*gaf*gaq*gs2sctw2*kp1m1*Qf*Qq*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * ( 64.D0 )                                  
      elmat2 = elmat2 + p1p3*p2p4*p3p4*gvf*gvq*gs2sctw2*kp3m1*kp4m1*Qf* 
     & Qq*Qf2*e4*dg12*dz12m2*s12mmz2 * ( 128.D0 )                       
      elmat2 = elmat2 + p1p3*p2p4*p3p4*gaf*gaq*gs2sctw2*kp3m1*kp4m1*Qf* 
     & Qq*Qf2*e4*dg12*dz12m2*s12mmz2 * (  - 128.D0 )                    
      elmat2 = elmat2 + p1p3*p2p3*gvf*gvq*gs2sctw2*kp3m1*Qf*Qq*Qf2*e4*  
     & dg12*dz12m2*s12mmz2 * (  - 128.D0 )                              
      elmat2 = elmat2 + p1p3*p2p3*gvf*gvq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp4 )                           
      elmat2 = elmat2 + p1p3*p2p3*gvf*gvq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp4 )                           
      elmat2 = elmat2 + p1p3*p2p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp4 )                        
      elmat2 = elmat2 + p1p3*p2p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp4 )                        
      elmat2 = elmat2 + p1p3*p2p3*gaf*gaq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp4 )                        
      elmat2 = elmat2 + p1p3*p2p3*gaf*gaq*gs2sctw2*kp2m1*kp3m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp4 )                        
      elmat2 = elmat2 + p1p3*p2p3*gaf*gaq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp4 )                        
      elmat2 = elmat2 + p1p3*p2p3*gaf*gaq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp4 )                        
      elmat2 = elmat2 + p1p3*p2p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp3m1*Qf2 
     & *Qq2*e4*dg34*dz12m2*s12mmz2 * (  - 64.D0 )                       
      elmat2 = elmat2 + p1p3*p2p3*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp3m1*Qf2 
     & *Qq2*e4*dg12*dz34m2*s34mmz2 * (  - 64.D0 )                       
      elmat2 = elmat2 + p1p3*p2p3*p2p4*gaf*gaq*gs2sctw2*kp2m1*kp3m1*Qf2 
     & *Qq2*e4*dg34*dz12m2*s12mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p3*p2p3*p2p4*gaf*gaq*gs2sctw2*kp2m1*kp3m1*Qf2 
     & *Qq2*e4*dg12*dz34m2*s34mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p3*p1p4*gvf*gvq*gs2sctw2*kp1m1*Qf*Qq*Qq2*e4*  
     & dg34*dz34m2*s34mmz2 * ( 128.D0 )                                 
      elmat2 = elmat2 + p1p3*p1p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * (  - 32.D0*kp2 )                        
      elmat2 = elmat2 + p1p3*p1p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * (  - 32.D0*kp2 )                        
      elmat2 = elmat2 + p1p3*p1p4*gvf*gvq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp2 )                           
      elmat2 = elmat2 + p1p3*p1p4*gvf*gvq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp2 )                           
      elmat2 = elmat2 + p1p3*p1p4*gaf*gaq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp2 )                           
      elmat2 = elmat2 + p1p3*p1p4*gaf*gaq*gs2sctw2*kp1m1*kp4m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp2 )                           
      elmat2 = elmat2 + p1p3*p1p4*gaf*gaq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2* 
     & e4*dg34*dz12m2*s12mmz2 * ( 32.D0*kp2 )                           
      elmat2 = elmat2 + p1p3*p1p4*gaf*gaq*gs2sctw2*kp1m1*kp3m1*Qf2*Qq2* 
     & e4*dg12*dz34m2*s34mmz2 * ( 32.D0*kp2 )                           
      elmat2 = elmat2 + p1p3*p1p4*p2p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*Qf2 
     & *Qq2*e4*dg34*dz12m2*s12mmz2 * (  - 64.D0 )                       
      elmat2 = elmat2 + p1p3*p1p4*p2p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*Qf2 
     & *Qq2*e4*dg12*dz34m2*s34mmz2 * (  - 64.D0 )                       
      elmat2 = elmat2 + p1p3*p1p4*p2p4*gaf*gaq*gs2sctw2*kp1m1*kp4m1*Qf2 
     & *Qq2*e4*dg34*dz12m2*s12mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p3*p1p4*p2p4*gaf*gaq*gs2sctw2*kp1m1*kp4m1*Qf2 
     & *Qq2*e4*dg12*dz34m2*s34mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p3*p1p4*p2p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*Qf2 
     & *Qq2*e4*dg34*dz12m2*s12mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p3*p1p4*p2p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*Qf2 
     & *Qq2*e4*dg12*dz34m2*s34mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p3*p1p4*p2p3*gaf*gaq*gs2sctw2*kp1m1*kp3m1*Qf2 
     & *Qq2*e4*dg34*dz12m2*s12mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p3*p1p4*p2p3*gaf*gaq*gs2sctw2*kp1m1*kp3m1*Qf2 
     & *Qq2*e4*dg12*dz34m2*s34mmz2 * ( 64.D0 )                          
      elmat2 = elmat2 + p1p2*gvf*gvq*gs2sctw2*Qf*Qq*mf4*Qf2*e4*dg12*    
     & dz12m2*s12mmz2 * (  - 64.D0*kp4m2 - 64.D0*kp3m2 )                
      elmat2 = elmat2 + p1p2*gvf*gvq*gs2sctw2*Qf*Qq*mq2*mf2*Qq2*e4*dg34 
     & *dz34m2*s34mmz2 * (  - 64.D0*kp2m2 - 64.D0*kp1m2 )               
      elmat2 = elmat2 + p1p2*gvf*gvq*gs2sctw2*kp2m1*Qf*Qq*mf2*Qq2*e4*   
     & dg34*dz34m2*s34mmz2 * (  - 128.D0 )                              
      elmat2 = elmat2 + p1p2*gvf*gvq*gs2sctw2*kp1m1*Qf*Qq*mf2*Qq2*e4*   
     & dg34*dz34m2*s34mmz2 * (  - 128.D0 )                              
      elmat2 = elmat2 + p1p2*gvf*gvq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*mq2*mf2 
     & *Qq2*e4*dg34*dz34m2*s34mmz2 * ( 256.D0 )                         
      elmat2 = elmat2 + p1p2*gaf*gaq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + p1p2*gaf*gaq*gs2sctw2*kp4m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + p1p2*gaf*gaq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg34 
     & *dz12m2*s12mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + p1p2*gaf*gaq*gs2sctw2*kp3m1*mf2*Qf2*Qq2*e4*dg12 
     & *dz34m2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + p1p2*p3p4*gvf*gvq*gs2sctw2*kp3m1*kp4m1*Qf*Qq*   
     & mf2*Qf2*e4*dg12*dz12m2*s12mmz2 * ( 128.D0 )                      
      elmat2 = elmat2 + p1p2*p3p4*gvf*gvq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*   
     & mq2*Qq2*e4*dg34*dz34m2*s34mmz2 * ( 128.D0 )                      
      elmat2 = elmat2 + p1p2*p2p4*gvf*gvq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p2*p2p4*gvf*gvq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p2*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*mf2*Qf2* 
     & Qq2*e4*dg34*dz12m2*s12mmz2 * ( 64.D0 )                           
      elmat2 = elmat2 + p1p2*p2p4*gvf*gvq*gs2sctw2*kp2m1*kp4m1*mf2*Qf2* 
     & Qq2*e4*dg12*dz34m2*s34mmz2 * ( 64.D0 )                           
      elmat2 = elmat2 + p1p2*p2p4*gvf*gvq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*   
     & Qq2*e4*dg34*dz34m2*s34mmz2 * (  - 64.D0*kp3 )                    
      elmat2 = elmat2 + p1p2*p2p4*gaf*gaq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p2p4*gaf*gaq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p2p4*gaf*gaq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*   
     & Qq2*e4*dg34*dz34m2*s34mmz2 * ( 64.D0*kp3 )                       
      elmat2 = elmat2 + p1p2*p2p3*gvf*gvq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p2p3*gvf*gvq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p2p3*gvf*gvq*gs2sctw2*kp2m1*kp3m1*mf2*Qf2* 
     & Qq2*e4*dg34*dz12m2*s12mmz2 * (  - 64.D0 )                        
      elmat2 = elmat2 + p1p2*p2p3*gvf*gvq*gs2sctw2*kp2m1*kp3m1*mf2*Qf2* 
     & Qq2*e4*dg12*dz34m2*s34mmz2 * (  - 64.D0 )                        
      elmat2 = elmat2 + p1p2*p2p3*gvf*gvq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*   
     & Qq2*e4*dg34*dz34m2*s34mmz2 * (  - 64.D0*kp4 )                    
      elmat2 = elmat2 + p1p2*p2p3*gaf*gaq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p2p3*gaf*gaq*gs2sctw2*kp2m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p2p3*gaf*gaq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*   
     & Qq2*e4*dg34*dz34m2*s34mmz2 * (  - 64.D0*kp4 )                    
      elmat2 = elmat2 + p1p2*p1p4*gvf*gvq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p1p4*gvf*gvq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p1p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*mf2*Qf2* 
     & Qq2*e4*dg34*dz12m2*s12mmz2 * (  - 64.D0 )                        
      elmat2 = elmat2 + p1p2*p1p4*gvf*gvq*gs2sctw2*kp1m1*kp4m1*mf2*Qf2* 
     & Qq2*e4*dg12*dz34m2*s34mmz2 * (  - 64.D0 )                        
      elmat2 = elmat2 + p1p2*p1p4*gvf*gvq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*   
     & Qq2*e4*dg34*dz34m2*s34mmz2 * (  - 64.D0*kp3 )                    
      elmat2 = elmat2 + p1p2*p1p4*gaf*gaq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p1p4*gaf*gaq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p1p4*gaf*gaq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*   
     & Qq2*e4*dg34*dz34m2*s34mmz2 * (  - 64.D0*kp3 )                    
      elmat2 = elmat2 + p1p2*p1p4*p2p3*gvf*gvq*gs2sctw2*kp1m1*kp2m1*Qf* 
     & Qq*Qq2*e4*dg34*dz34m2*s34mmz2 * ( 128.D0 )                       
      elmat2 = elmat2 + p1p2*p1p4*p2p3*gaf*gaq*gs2sctw2*kp1m1*kp2m1*Qf* 
     & Qq*Qq2*e4*dg34*dz34m2*s34mmz2 * ( 128.D0 )                       
      elmat2 = elmat2 + p1p2*p1p3*gvf*gvq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p2*p1p3*gvf*gvq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + p1p2*p1p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*mf2*Qf2* 
     & Qq2*e4*dg34*dz12m2*s12mmz2 * ( 64.D0 )                           
      elmat2 = elmat2 + p1p2*p1p3*gvf*gvq*gs2sctw2*kp1m1*kp3m1*mf2*Qf2* 
     & Qq2*e4*dg12*dz34m2*s34mmz2 * ( 64.D0 )                           
      elmat2 = elmat2 + p1p2*p1p3*gvf*gvq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*   
     & Qq2*e4*dg34*dz34m2*s34mmz2 * (  - 64.D0*kp4 )                    
      elmat2 = elmat2 + p1p2*p1p3*gaf*gaq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg34*dz12m2*s12mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p1p3*gaf*gaq*gs2sctw2*kp1m1*Qf2*Qq2*e4*    
     & dg12*dz34m2*s34mmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + p1p2*p1p3*gaf*gaq*gs2sctw2*kp1m1*kp2m1*Qf*Qq*   
     & Qq2*e4*dg34*dz34m2*s34mmz2 * ( 64.D0*kp4 )                       
      elmat2 = elmat2 + p1p2*p1p3*p2p4*gvf*gvq*gs2sctw2*kp1m1*kp2m1*Qf* 
     & Qq*Qq2*e4*dg34*dz34m2*s34mmz2 * ( 128.D0 )                       
      elmat2 = elmat2 + p1p2*p1p3*p2p4*gaf*gaq*gs2sctw2*kp1m1*kp2m1*Qf* 
     & Qq*Qq2*e4*dg34*dz34m2*s34mmz2 * (  - 128.D0 )                    
      elmat2 = elmat2 + e2*mq4*mf2*Qf2*Qq4*e4*dg342 * (  - 64.D0*kp2m2  
     &     - 64.D0*kp1m2 )                                              
      elmat2 = elmat2 + e2*mq2*mf4*Qf4*Qq2*e4*dg122 * (  - 64.D0*kp4m2  
     &     - 64.D0*kp3m2 )                                              
      elmat2 = elmat2 + e2*mq2*mf2*Qf4*Qq2*e4*dg122 * (  - 32.D0*kp4*   
     &    kp3m2 - 32.D0*kp3*kp4m2 )                                     
      elmat2 = elmat2 + e2*mq2*mf2*Qf2*Qq4*e4*dg342 * ( 32.D0*kp2*kp1m2 
     &     + 32.D0*kp1*kp2m2 )                                          
      elmat2 = elmat2 + e2*kp4m1*mq2*Qf4*Qq2*e4*dg122 * ( 32.D0*kp3 )   
      elmat2 = elmat2 + e2*kp4m1*mq2*mf2*Qf4*Qq2*e4*dg122 * (  - 32.D0  
     &     )                                                            
      elmat2 = elmat2 + e2*kp4m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*dg34 * ( 32.D
     &    0*kp2 - 32.D0*kp1 )                                           
      elmat2 = elmat2 + e2*kp3m1*mq2*Qf4*Qq2*e4*dg122 * ( 32.D0*kp4 )   
      elmat2 = elmat2 + e2*kp3m1*mq2*mf2*Qf4*Qq2*e4*dg122 * (  - 32.D0  
     &     )                                                            
      elmat2 = elmat2 + e2*kp3m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*dg34 * (  -  
     &    32.D0*kp2 + 32.D0*kp1 )                                       
      elmat2 = elmat2 + e2*kp3m1*kp4m1*mf2*Qf4*Qq2*e4*dg122 * (  - 64.D0
     &    *kp1*kp2 )                                                    
      elmat2 = elmat2 + e2*kp3m1*kp4m1*mq2*Qf4*Qq2*e4*dg122 * ( 64.D0*  
     &    p3p42 )                                                       
      elmat2 = elmat2 + e2*kp2m1*mf2*Qf2*Qq4*e4*dg342 * ( 32.D0*kp1 )   
      elmat2 = elmat2 + e2*kp2m1*mq2*mf2*Qf2*Qq4*e4*dg342 * ( 32.D0 )   
      elmat2 = elmat2 + e2*kp2m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*dg34 * ( 32.D
     &    0*kp4 - 32.D0*kp3 )                                           
      elmat2 = elmat2 + e2*kp1m1*mf2*Qf2*Qq4*e4*dg342 * ( 32.D0*kp2 )   
      elmat2 = elmat2 + e2*kp1m1*mq2*mf2*Qf2*Qq4*e4*dg342 * ( 32.D0 )   
      elmat2 = elmat2 + e2*kp1m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*dg34 * (  -  
     &    32.D0*kp4 + 32.D0*kp3 )                                       
      elmat2 = elmat2 + e2*kp1m1*kp2m1*mf2*Qf2*Qq4*e4*dg342 * ( 64.D0*  
     &    p1p22 )                                                       
      elmat2 = elmat2 + e2*kp1m1*kp2m1*mq2*Qf2*Qq4*e4*dg342 * (  - 64.D0
     &    *kp3*kp4 )                                                    
      elmat2 = elmat2 + e2*gs2sctw4*mq4*mf2*Qq2*gvf2*gvq2*dz34m2 * (    
     &     - 64.D0*kp2m2 - 64.D0*kp1m2 )                                
      elmat2 = elmat2 + e2*gs2sctw4*mq4*mf2*Qq2*gaq2*gvf2*dz34m2 * ( 64.
     &    D0*kp2m2 + 64.D0*kp1m2 )                                      
      elmat2 = elmat2 + e2*gs2sctw4*mq4*mf2*Qq2*gaf2*gvq2*dz34m2 * ( 64.
     &    D0*kp2m2 + 64.D0*kp1m2 )                                      
      elmat2 = elmat2 + e2*gs2sctw4*mq4*mf2*Qq2*gaf2*gaq2*dz34m2 * (    
     &     - 64.D0*kp2m2 - 64.D0*kp1m2 )                                
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf4*Qf2*gvf2*gvq2*dz12m2 * (    
     &     - 64.D0*kp4m2 - 64.D0*kp3m2 )                                
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf4*Qf2*gaq2*gvf2*dz12m2 * ( 64.
     &    D0*kp4m2 + 64.D0*kp3m2 )                                      
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf4*Qf2*gaf2*gvq2*dz12m2 * ( 64.
     &    D0*kp4m2 + 64.D0*kp3m2 )                                      
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf4*Qf2*gaf2*gaq2*dz12m2 * (    
     &     - 64.D0*kp4m2 - 64.D0*kp3m2 )                                
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf2*Qq2*gvf2*gvq2*dz34m2 * ( 32.
     &    D0*kp2*kp1m2 + 32.D0*kp1*kp2m2 )                              
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf2*Qq2*gaq2*gvf2*dz34m2 * ( 32.
     &    D0*kp2*kp1m2 + 32.D0*kp1*kp2m2 )                              
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf2*Qq2*gaf2*gvq2*dz34m2 * (    
     &     - 32.D0*kp2*kp1m2 - 32.D0*kp1*kp2m2 )                        
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf2*Qq2*gaf2*gaq2*dz34m2 * (    
     &     - 32.D0*kp2*kp1m2 - 32.D0*kp1*kp2m2 )                        
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf2*Qf2*gvf2*gvq2*dz12m2 * (    
     &     - 32.D0*kp4*kp3m2 - 32.D0*kp3*kp4m2 )                        
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf2*Qf2*gaq2*gvf2*dz12m2 * ( 32.
     &    D0*kp4*kp3m2 + 32.D0*kp3*kp4m2 )                              
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf2*Qf2*gaf2*gvq2*dz12m2 * (    
     &     - 32.D0*kp4*kp3m2 - 32.D0*kp3*kp4m2 )                        
      elmat2 = elmat2 + e2*gs2sctw4*mq2*mf2*Qf2*gaf2*gaq2*dz12m2 * ( 32.
     &    D0*kp4*kp3m2 + 32.D0*kp3*kp4m2 )                              
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*mq2*Qf2*gvf2*gvq2*dz12m2 * (  
     &    32.D0*kp3 )                                                   
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*mq2*Qf2*gaq2*gvf2*dz12m2 * (  
     &     - 32.D0*kp3 )                                                
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*mq2*Qf2*gaf2*gvq2*dz12m2 * (  
     &    32.D0*kp3 )                                                   
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*mq2*Qf2*gaf2*gaq2*dz12m2 * (  
     &     - 32.D0*kp3 )                                                
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*mq2*mf2*Qf2*gvf2*gvq2*dz12m2  
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*mq2*mf2*Qf2*gaq2*gvf2*dz12m2  
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*mq2*mf2*Qf2*gaf2*gvq2*dz12m2  
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*mq2*mf2*Qf2*gaf2*gaq2*dz12m2  
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*Qf*Qq*mf2*gvf2*gvq2*dz12m2*   
     & dz34m2*gzmz2 * ( 32.D0*kp2 - 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*Qf*Qq*mf2*gvf2*gvq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 - 32.D0*kp1 )               
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*Qf*Qq*mf2*gaq2*gvf2*dz12m2*   
     & dz34m2*gzmz2 * ( 32.D0*kp2 - 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*Qf*Qq*mf2*gaq2*gvf2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 - 32.D0*kp1 )               
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gvq2*dz12m2*   
     & dz34m2*gzmz2 * (  - 32.D0*kp2 + 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gvq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 + 32.D0*kp1 )            
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gaq2*dz12m2*   
     & dz34m2*gzmz2 * (  - 32.D0*kp2 + 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gaq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 + 32.D0*kp1 )            
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*mq2*Qf2*gvf2*gvq2*dz12m2 * (  
     &    32.D0*kp4 )                                                   
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*mq2*Qf2*gaq2*gvf2*dz12m2 * (  
     &     - 32.D0*kp4 )                                                
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*mq2*Qf2*gaf2*gvq2*dz12m2 * (  
     &    32.D0*kp4 )                                                   
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*mq2*Qf2*gaf2*gaq2*dz12m2 * (  
     &     - 32.D0*kp4 )                                                
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*mq2*mf2*Qf2*gvf2*gvq2*dz12m2  
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*mq2*mf2*Qf2*gaq2*gvf2*dz12m2  
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*mq2*mf2*Qf2*gaf2*gvq2*dz12m2  
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*mq2*mf2*Qf2*gaf2*gaq2*dz12m2  
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*Qf*Qq*mf2*gvf2*gvq2*dz12m2*   
     & dz34m2*gzmz2 * (  - 32.D0*kp2 + 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*Qf*Qq*mf2*gvf2*gvq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 + 32.D0*kp1 )            
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*Qf*Qq*mf2*gaq2*gvf2*dz12m2*   
     & dz34m2*gzmz2 * (  - 32.D0*kp2 + 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*Qf*Qq*mf2*gaq2*gvf2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 + 32.D0*kp1 )            
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gvq2*dz12m2*   
     & dz34m2*gzmz2 * ( 32.D0*kp2 - 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gvq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 - 32.D0*kp1 )               
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gaq2*dz12m2*   
     & dz34m2*gzmz2 * ( 32.D0*kp2 - 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gaq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 - 32.D0*kp1 )               
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*kp4m1*mf2*Qf2*gvf2*gvq2*      
     & dz12m2 * (  - 64.D0*kp1*kp2 )                                    
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*kp4m1*mf2*Qf2*gaq2*gvf2*      
     & dz12m2 * (  - 64.D0*kp1*kp2 )                                    
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*kp4m1*mf2*Qf2*gaf2*gvq2*      
     & dz12m2 * ( 64.D0*kp1*kp2 )                                       
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*kp4m1*mf2*Qf2*gaf2*gaq2*      
     & dz12m2 * ( 64.D0*kp1*kp2 )                                       
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*kp4m1*mq2*Qf2*gvf2*gvq2*      
     & dz12m2 * ( 64.D0*p3p42 )                                         
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*kp4m1*mq2*Qf2*gaq2*gvf2*      
     & dz12m2 * (  - 64.D0*p3p42 )                                      
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*kp4m1*mq2*Qf2*gaf2*gvq2*      
     & dz12m2 * ( 64.D0*p3p42 )                                         
      elmat2 = elmat2 + e2*gs2sctw4*kp3m1*kp4m1*mq2*Qf2*gaf2*gaq2*      
     & dz12m2 * (  - 64.D0*p3p42 )                                      
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*mf2*Qq2*gvf2*gvq2*dz34m2 * (  
     &    32.D0*kp1 )                                                   
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*mf2*Qq2*gaq2*gvf2*dz34m2 * (  
     &    32.D0*kp1 )                                                   
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*mf2*Qq2*gaf2*gvq2*dz34m2 * (  
     &     - 32.D0*kp1 )                                                
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*mf2*Qq2*gaf2*gaq2*dz34m2 * (  
     &     - 32.D0*kp1 )                                                
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*mq2*mf2*Qq2*gvf2*gvq2*dz34m2  
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*mq2*mf2*Qq2*gaq2*gvf2*dz34m2  
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*mq2*mf2*Qq2*gaf2*gvq2*dz34m2  
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*mq2*mf2*Qq2*gaf2*gaq2*dz34m2  
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*Qf*Qq*mq2*gvf2*gvq2*dz12m2*   
     & dz34m2*gzmz2 * ( 32.D0*kp4 - 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*Qf*Qq*mq2*gvf2*gvq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 - 32.D0*kp3 )               
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*Qf*Qq*mq2*gaq2*gvf2*dz12m2*   
     & dz34m2*gzmz2 * (  - 32.D0*kp4 + 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*Qf*Qq*mq2*gaq2*gvf2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 + 32.D0*kp3 )            
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gvq2*dz12m2*   
     & dz34m2*gzmz2 * ( 32.D0*kp4 - 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gvq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 - 32.D0*kp3 )               
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gaq2*dz12m2*   
     & dz34m2*gzmz2 * (  - 32.D0*kp4 + 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gaq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 + 32.D0*kp3 )            
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*mf2*Qq2*gvf2*gvq2*dz34m2 * (  
     &    32.D0*kp2 )                                                   
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*mf2*Qq2*gaq2*gvf2*dz34m2 * (  
     &    32.D0*kp2 )                                                   
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*mf2*Qq2*gaf2*gvq2*dz34m2 * (  
     &     - 32.D0*kp2 )                                                
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*mf2*Qq2*gaf2*gaq2*dz34m2 * (  
     &     - 32.D0*kp2 )                                                
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*mq2*mf2*Qq2*gvf2*gvq2*dz34m2  
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*mq2*mf2*Qq2*gaq2*gvf2*dz34m2  
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*mq2*mf2*Qq2*gaf2*gvq2*dz34m2  
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*mq2*mf2*Qq2*gaf2*gaq2*dz34m2  
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*Qf*Qq*mq2*gvf2*gvq2*dz12m2*   
     & dz34m2*gzmz2 * (  - 32.D0*kp4 + 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*Qf*Qq*mq2*gvf2*gvq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 + 32.D0*kp3 )            
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*Qf*Qq*mq2*gaq2*gvf2*dz12m2*   
     & dz34m2*gzmz2 * ( 32.D0*kp4 - 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*Qf*Qq*mq2*gaq2*gvf2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 - 32.D0*kp3 )               
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gvq2*dz12m2*   
     & dz34m2*gzmz2 * (  - 32.D0*kp4 + 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gvq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 + 32.D0*kp3 )            
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gaq2*dz12m2*   
     & dz34m2*gzmz2 * ( 32.D0*kp4 - 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gaq2*dz12m2*   
     & dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 - 32.D0*kp3 )               
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*kp2m1*mf2*Qq2*gvf2*gvq2*      
     & dz34m2 * ( 64.D0*p1p22 )                                         
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*kp2m1*mf2*Qq2*gaq2*gvf2*      
     & dz34m2 * ( 64.D0*p1p22 )                                         
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*kp2m1*mf2*Qq2*gaf2*gvq2*      
     & dz34m2 * (  - 64.D0*p1p22 )                                      
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*kp2m1*mf2*Qq2*gaf2*gaq2*      
     & dz34m2 * (  - 64.D0*p1p22 )                                      
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*kp2m1*mq2*Qq2*gvf2*gvq2*      
     & dz34m2 * (  - 64.D0*kp3*kp4 )                                    
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*kp2m1*mq2*Qq2*gaq2*gvf2*      
     & dz34m2 * ( 64.D0*kp3*kp4 )                                       
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*kp2m1*mq2*Qq2*gaf2*gvq2*      
     & dz34m2 * (  - 64.D0*kp3*kp4 )                                    
      elmat2 = elmat2 + e2*gs2sctw4*kp1m1*kp2m1*mq2*Qq2*gaf2*gaq2*      
     & dz34m2 * ( 64.D0*kp3*kp4 )                                       
      elmat2 = elmat2 + e2*p3p4*mq4*Qf2*Qq4*e4*dg342 * (  - 32.D0*kp2m2 
     &     - 32.D0*kp1m2 )                                              
      elmat2 = elmat2 + e2*p3p4*mq2*mf2*Qf4*Qq2*e4*dg122 * (  - 32.D0*  
     &    kp4m2 - 32.D0*kp3m2 )                                         
      elmat2 = elmat2 + e2*p3p4*kp4m1*mq2*Qf4*Qq2*e4*dg122 * ( 64.D0 )  
      elmat2 = elmat2 + e2*p3p4*kp3m1*mq2*Qf4*Qq2*e4*dg122 * ( 64.D0 )  
      elmat2 = elmat2 + e2*p3p4*kp3m1*kp4m1*mq2*mf2*Qf4*Qq2*e4*dg122    
     &  * ( 128.D0 )                                                    
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*mq4*Qq2*gvf2*gvq2*dz34m2 * (   
     &     - 32.D0*kp2m2 - 32.D0*kp1m2 )                                
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*mq4*Qq2*gaq2*gvf2*dz34m2 * (   
     &    32.D0*kp2m2 + 32.D0*kp1m2 )                                   
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*mq4*Qq2*gaf2*gvq2*dz34m2 * (   
     &     - 32.D0*kp2m2 - 32.D0*kp1m2 )                                
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*mq4*Qq2*gaf2*gaq2*dz34m2 * (   
     &    32.D0*kp2m2 + 32.D0*kp1m2 )                                   
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*mq2*mf2*Qf2*gvf2*gvq2*dz12m2   
     &  * (  - 32.D0*kp4m2 - 32.D0*kp3m2 )                              
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*mq2*mf2*Qf2*gaq2*gvf2*dz12m2   
     &  * ( 32.D0*kp4m2 + 32.D0*kp3m2 )                                 
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*mq2*mf2*Qf2*gaf2*gvq2*dz12m2   
     &  * (  - 32.D0*kp4m2 - 32.D0*kp3m2 )                              
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*mq2*mf2*Qf2*gaf2*gaq2*dz12m2   
     &  * ( 32.D0*kp4m2 + 32.D0*kp3m2 )                                 
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp4m1*mq2*Qf2*gvf2*gvq2*dz12m2 
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp4m1*mq2*Qf2*gaq2*gvf2*dz12m2 
     &  * (  - 64.D0 )                                                  
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp4m1*mq2*Qf2*gaf2*gvq2*dz12m2 
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp4m1*mq2*Qf2*gaf2*gaq2*dz12m2 
     &  * (  - 64.D0 )                                                  
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp3m1*mq2*Qf2*gvf2*gvq2*dz12m2 
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp3m1*mq2*Qf2*gaq2*gvf2*dz12m2 
     &  * (  - 64.D0 )                                                  
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp3m1*mq2*Qf2*gaf2*gvq2*dz12m2 
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp3m1*mq2*Qf2*gaf2*gaq2*dz12m2 
     &  * (  - 64.D0 )                                                  
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp3m1*kp4m1*mq2*mf2*Qf2*gvf2*  
     & gvq2*dz12m2 * ( 128.D0 )                                         
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp3m1*kp4m1*mq2*mf2*Qf2*gaq2*  
     & gvf2*dz12m2 * (  - 128.D0 )                                      
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp3m1*kp4m1*mq2*mf2*Qf2*gaf2*  
     & gvq2*dz12m2 * (  - 128.D0 )                                      
      elmat2 = elmat2 + e2*p3p4*gs2sctw4*kp3m1*kp4m1*mq2*mf2*Qf2*gaf2*  
     & gaq2*dz12m2 * ( 128.D0 )                                         
      elmat2 = elmat2 + e2*p2p4*mf2*Qf4*Qq2*e4*dg122 * (  - 32.D0*kp1*  
     &    kp3m2 )                                                       
      elmat2 = elmat2 + e2*p2p4*mq2*Qf2*Qq4*e4*dg342 * ( 32.D0*kp3*     
     &    kp1m2 )                                                       
      elmat2 = elmat2 + e2*p2p4*Qf*Qq*Qf2*Qq2*e4*dg12*dg34 * (  - 64.D0 
     &     )                                                            
      elmat2 = elmat2 + e2*p2p4*kp4m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*dg34    
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p2p4*kp3m1*Qf4*Qq2*e4*dg122 * ( 32.D0*kp1 )  
      elmat2 = elmat2 + e2*p2p4*kp3m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*dg34    
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p2p4*kp2m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*dg34    
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p2p4*kp2m1*kp4m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*  
     & dg34 * (  - 32.D0*kp1 )                                          
      elmat2 = elmat2 + e2*p2p4*kp2m1*kp4m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*  
     & dg34 * ( 32.D0*kp3 )                                             
      elmat2 = elmat2 + e2*p2p4*kp2m1*kp4m1*Qf*Qq*mq2*mf2*Qf2*Qq2*e4*   
     & dg12*dg34 * ( 128.D0 )                                           
      elmat2 = elmat2 + e2*p2p4*kp1m1*Qf2*Qq4*e4*dg342 * ( 32.D0*kp3 )  
      elmat2 = elmat2 + e2*p2p4*kp1m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*dg34    
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p2p4*kp1m1*kp3m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34  
     &  * ( 64.D0*p1p32 )                                               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*mf2*Qf2*gvf2*gvq2*dz12m2 * (   
     &     - 32.D0*kp1*kp3m2 )                                          
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*mf2*Qf2*gaq2*gvf2*dz12m2 * (   
     &     - 32.D0*kp1*kp3m2 )                                          
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*mf2*Qf2*gaf2*gvq2*dz12m2 * (   
     &     - 32.D0*kp1*kp3m2 )                                          
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*mf2*Qf2*gaf2*gaq2*dz12m2 * (   
     &     - 32.D0*kp1*kp3m2 )                                          
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*mq2*Qq2*gvf2*gvq2*dz34m2 * (   
     &    32.D0*kp3*kp1m2 )                                             
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*mq2*Qq2*gaq2*gvf2*dz34m2 * (   
     &    32.D0*kp3*kp1m2 )                                             
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*mq2*Qq2*gaf2*gvq2*dz34m2 * (   
     &    32.D0*kp3*kp1m2 )                                             
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*mq2*Qq2*gaf2*gaq2*dz34m2 * (   
     &    32.D0*kp3*kp1m2 )                                             
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*Qf*Qq*gvf2*gvq2*dz12m2*dz34m2* 
     & gzmz2 * (  - 64.D0 )                                             
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*Qf*Qq*gvf2*gvq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*Qf*Qq*gaq2*gvf2*dz12m2*dz34m2* 
     & gzmz2 * (  - 64.D0 )                                             
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*Qf*Qq*gaq2*gvf2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*Qf*Qq*gaf2*gvq2*dz12m2*dz34m2* 
     & gzmz2 * (  - 64.D0 )                                             
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*Qf*Qq*gaf2*gvq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*Qf*Qq*gaf2*gaq2*dz12m2*dz34m2* 
     & gzmz2 * (  - 64.D0 )                                             
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*Qf*Qq*gaf2*gaq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf2*gvf2*gvq2*dz12m2     
     &  * ( 32.D0*kp1 )                                                 
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf2*gaq2*gvf2*dz12m2     
     &  * ( 32.D0*kp1 )                                                 
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf2*gaf2*gvq2*dz12m2     
     &  * ( 32.D0*kp1 )                                                 
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf2*gaf2*gaq2*dz12m2     
     &  * ( 32.D0*kp1 )                                                 
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*gvf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*gvf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp1 )            
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*gaq2*    
     & gvf2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*gaq2*    
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp1 )            
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*gaf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*gaf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp1 )               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*gaf2*    
     & gaq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*gaf2*    
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp1 )               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*gvf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*gvf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp3 )               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*gaq2*    
     & gvf2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*gaq2*    
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp3 )            
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*gaf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*gaf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp3 )               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*gaf2*    
     & gaq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*gaf2*    
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp3 )            
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*mf2*gvf2 
     & *gvq2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                           
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*mf2*gvf2 
     & *gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                 
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*mf2*gaq2 
     & *gvf2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                        
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*mf2*gaq2 
     & *gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )              
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*mf2*gaf2 
     & *gvq2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                        
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*mf2*gaf2 
     & *gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )              
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*mf2*gaf2 
     & *gaq2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                           
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*mf2*gaf2 
     & *gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                 
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qq2*gvf2*gvq2*dz34m2     
     &  * ( 32.D0*kp3 )                                                 
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qq2*gaq2*gvf2*dz34m2     
     &  * ( 32.D0*kp3 )                                                 
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qq2*gaf2*gvq2*dz34m2     
     &  * ( 32.D0*kp3 )                                                 
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qq2*gaf2*gaq2*dz34m2     
     &  * ( 32.D0*kp3 )                                                 
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gvf2*gvq2*   
     & dz12m2*dz34m2*gzmz2 * ( 64.D0*p1p32 )                            
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gvf2*gvq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0*p1p32 )                  
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaq2*gvf2*   
     & dz12m2*dz34m2*gzmz2 * ( 64.D0*p1p32 )                            
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaq2*gvf2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0*p1p32 )                  
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*gvq2*   
     & dz12m2*dz34m2*gzmz2 * ( 64.D0*p1p32 )                            
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*gvq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0*p1p32 )                  
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*gaq2*   
     & dz12m2*dz34m2*gzmz2 * ( 64.D0*p1p32 )                            
      elmat2 = elmat2 + e2*p2p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*gaq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0*p1p32 )                  
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*mf2*Qf2*dz12m2 
     &  * ( 128.D0*kp1*kp3m2 )                                          
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*mq2*Qq2*dz34m2 
     &  * (  - 128.D0*kp3*kp1m2 )                                       
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*Qf*Qq*dz12m2*  
     & dz34m2*gzmz2 * ( 256.D0 )                                        
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*Qf*Qq*dz12m2*  
     & dz34m2*s12mmz2*s34mmz2 * ( 256.D0 )                              
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf2*     
     & dz12m2 * (  - 128.D0*kp1 )                                       
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf*Qq*   
     & mf2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                          
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf*Qq*   
     & mf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qq2*     
     & dz34m2 * (  - 128.D0*kp3 )                                       
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf*Qq*   
     & mq2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                             
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf*Qq*   
     & mq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                   
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*kp3m1*Qf 
     & *Qq*dz12m2*dz34m2*gzmz2 * (  - 256.D0*p1p32 )                    
      elmat2 = elmat2 + e2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*kp3m1*Qf 
     & *Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 256.D0*p1p32 )          
      elmat2 = elmat2 + e2*p2p4*p3p4*kp4m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p2p4*p3p4*kp3m1*kp4m1*Qf4*Qq2*e4*dg122 * (   
     &    32.D0*kp1 )                                                   
      elmat2 = elmat2 + e2*p2p4*p3p4*kp2m1*kp4m1*Qf*Qq*mq2*Qf2*Qq2*e4*  
     & dg12*dg34 * ( 64.D0 )                                            
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gvf2*gvq2 
     & *dz12m2 * ( 32.D0*kp1 )                                          
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaq2*gvf2 
     & *dz12m2 * ( 32.D0*kp1 )                                          
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2*gvq2 
     & *dz12m2 * ( 32.D0*kp1 )                                          
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2*gaq2 
     & *dz12m2 * ( 32.D0*kp1 )                                          
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*    
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*    
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*    
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*    
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*    
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*    
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*    
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p2p4*p3p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mq2*    
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p2p4*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                              
      elmat2 = elmat2 + e2*p2p4*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                    
      elmat2 = elmat2 + e2*p2p4*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*    
     & kp4m1*Qf2*dz12m2 * (  - 128.D0*kp1 )                             
      elmat2 = elmat2 + e2*p2p3*mf2*Qf4*Qq2*e4*dg122 * (  - 32.D0*kp1*  
     &    kp4m2 )                                                       
      elmat2 = elmat2 + e2*p2p3*mq2*Qf2*Qq4*e4*dg342 * ( 32.D0*kp4*     
     &    kp1m2 )                                                       
      elmat2 = elmat2 + e2*p2p3*Qf*Qq*Qf2*Qq2*e4*dg12*dg34 * ( 64.D0 )  
      elmat2 = elmat2 + e2*p2p3*kp4m1*Qf4*Qq2*e4*dg122 * ( 32.D0*kp1 )  
      elmat2 = elmat2 + e2*p2p3*kp4m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*dg34    
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p2p3*kp3m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*dg34    
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p2p3*kp2m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*dg34    
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p2p3*kp2m1*kp3m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*  
     & dg34 * ( 32.D0*kp1 )                                             
      elmat2 = elmat2 + e2*p2p3*kp2m1*kp3m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*  
     & dg34 * (  - 32.D0*kp4 )                                          
      elmat2 = elmat2 + e2*p2p3*kp2m1*kp3m1*Qf*Qq*mq2*mf2*Qf2*Qq2*e4*   
     & dg12*dg34 * (  - 128.D0 )                                        
      elmat2 = elmat2 + e2*p2p3*kp1m1*Qf2*Qq4*e4*dg342 * ( 32.D0*kp4 )  
      elmat2 = elmat2 + e2*p2p3*kp1m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*dg34    
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p2p3*kp1m1*kp4m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34  
     &  * (  - 64.D0*p1p42 )                                            
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*mf2*Qf2*gvf2*gvq2*dz12m2 * (   
     &     - 32.D0*kp1*kp4m2 )                                          
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*mf2*Qf2*gaq2*gvf2*dz12m2 * (   
     &     - 32.D0*kp1*kp4m2 )                                          
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*mf2*Qf2*gaf2*gvq2*dz12m2 * (   
     &     - 32.D0*kp1*kp4m2 )                                          
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*mf2*Qf2*gaf2*gaq2*dz12m2 * (   
     &     - 32.D0*kp1*kp4m2 )                                          
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*mq2*Qq2*gvf2*gvq2*dz34m2 * (   
     &    32.D0*kp4*kp1m2 )                                             
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*mq2*Qq2*gaq2*gvf2*dz34m2 * (   
     &    32.D0*kp4*kp1m2 )                                             
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*mq2*Qq2*gaf2*gvq2*dz34m2 * (   
     &    32.D0*kp4*kp1m2 )                                             
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*mq2*Qq2*gaf2*gaq2*dz34m2 * (   
     &    32.D0*kp4*kp1m2 )                                             
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*Qf*Qq*gvf2*gvq2*dz12m2*dz34m2* 
     & gzmz2 * ( 64.D0 )                                                
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*Qf*Qq*gvf2*gvq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*Qf*Qq*gaq2*gvf2*dz12m2*dz34m2* 
     & gzmz2 * ( 64.D0 )                                                
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*Qf*Qq*gaq2*gvf2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*Qf*Qq*gaf2*gvq2*dz12m2*dz34m2* 
     & gzmz2 * ( 64.D0 )                                                
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*Qf*Qq*gaf2*gvq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*Qf*Qq*gaf2*gaq2*dz12m2*dz34m2* 
     & gzmz2 * ( 64.D0 )                                                
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*Qf*Qq*gaf2*gaq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf2*gvf2*gvq2*dz12m2     
     &  * ( 32.D0*kp1 )                                                 
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf2*gaq2*gvf2*dz12m2     
     &  * ( 32.D0*kp1 )                                                 
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf2*gaf2*gvq2*dz12m2     
     &  * ( 32.D0*kp1 )                                                 
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf2*gaf2*gaq2*dz12m2     
     &  * ( 32.D0*kp1 )                                                 
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*gvf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*gvf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp1 )               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*gaq2*    
     & gvf2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*gaq2*    
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp1 )               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*gaf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*gaf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp1 )            
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*gaf2*    
     & gaq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*gaf2*    
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp1 )            
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*gvf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp4 )                      
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*gvf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 )            
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*gaq2*    
     & gvf2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp4 )                         
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*gaq2*    
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 )               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*gaf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp4 )                      
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*gaf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 )            
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*gaf2*    
     & gaq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp4 )                         
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*gaf2*    
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 )               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*mf2*gvf2 
     & *gvq2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                        
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*mf2*gvf2 
     & *gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )              
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*mf2*gaq2 
     & *gvf2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                           
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*mf2*gaq2 
     & *gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                 
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*mf2*gaf2 
     & *gvq2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                           
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*mf2*gaf2 
     & *gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                 
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*mf2*gaf2 
     & *gaq2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                        
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*mf2*gaf2 
     & *gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )              
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qq2*gvf2*gvq2*dz34m2     
     &  * ( 32.D0*kp4 )                                                 
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qq2*gaq2*gvf2*dz34m2     
     &  * ( 32.D0*kp4 )                                                 
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qq2*gaf2*gvq2*dz34m2     
     &  * ( 32.D0*kp4 )                                                 
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qq2*gaf2*gaq2*dz34m2     
     &  * ( 32.D0*kp4 )                                                 
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gvf2*gvq2*   
     & dz12m2*dz34m2*gzmz2 * (  - 64.D0*p1p42 )                         
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gvf2*gvq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0*p1p42 )               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaq2*gvf2*   
     & dz12m2*dz34m2*gzmz2 * (  - 64.D0*p1p42 )                         
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaq2*gvf2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0*p1p42 )               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*gvq2*   
     & dz12m2*dz34m2*gzmz2 * (  - 64.D0*p1p42 )                         
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*gvq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0*p1p42 )               
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*gaq2*   
     & dz12m2*dz34m2*gzmz2 * (  - 64.D0*p1p42 )                         
      elmat2 = elmat2 + e2*p2p3*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*gaq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0*p1p42 )               
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*mf2*Qf2*dz12m2 
     &  * (  - 128.D0*kp1*kp4m2 )                                       
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*mq2*Qq2*dz34m2 
     &  * ( 128.D0*kp4*kp1m2 )                                          
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*Qf*Qq*dz12m2*  
     & dz34m2*gzmz2 * ( 256.D0 )                                        
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*Qf*Qq*dz12m2*  
     & dz34m2*s12mmz2*s34mmz2 * ( 256.D0 )                              
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf2*     
     & dz12m2 * ( 128.D0*kp1 )                                          
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf*Qq*   
     & mf2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                          
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf*Qq*   
     & mf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qq2*     
     & dz34m2 * ( 128.D0*kp4 )                                          
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf*Qq*   
     & mq2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                             
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf*Qq*   
     & mq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                   
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*kp4m1*Qf 
     & *Qq*dz12m2*dz34m2*gzmz2 * (  - 256.D0*p1p42 )                    
      elmat2 = elmat2 + e2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*kp4m1*Qf 
     & *Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 256.D0*p1p42 )          
      elmat2 = elmat2 + e2*p2p3*p3p4*kp3m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p2p3*p3p4*kp3m1*kp4m1*Qf4*Qq2*e4*dg122 * (   
     &    32.D0*kp1 )                                                   
      elmat2 = elmat2 + e2*p2p3*p3p4*kp2m1*kp3m1*Qf*Qq*mq2*Qf2*Qq2*e4*  
     & dg12*dg34 * (  - 64.D0 )                                         
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gvf2*gvq2 
     & *dz12m2 * ( 32.D0*kp1 )                                          
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaq2*gvf2 
     & *dz12m2 * ( 32.D0*kp1 )                                          
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2*gvq2 
     & *dz12m2 * ( 32.D0*kp1 )                                          
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2*gaq2 
     & *dz12m2 * ( 32.D0*kp1 )                                          
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*    
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*    
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*    
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*    
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*    
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*    
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*    
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p2p3*p3p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mq2*    
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p2p3*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                              
      elmat2 = elmat2 + e2*p2p3*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                    
      elmat2 = elmat2 + e2*p2p3*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*    
     & kp4m1*Qf2*dz12m2 * ( 128.D0*kp1 )                                
      elmat2 = elmat2 + e2*p2p3*p2p4*kp2m1*Qf2*Qq4*e4*dg342 * ( 64.D0 ) 
      elmat2 = elmat2 + e2*p2p3*p2p4*kp2m1*kp4m1*Qf*Qq*Qf2*Qq2*e4*dg12* 
     & dg34 * ( 32.D0*kp1 )                                             
      elmat2 = elmat2 + e2*p2p3*p2p4*kp2m1*kp3m1*Qf*Qq*Qf2*Qq2*e4*dg12* 
     & dg34 * (  - 32.D0*kp1 )                                          
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*Qq2*gvf2*gvq2*      
     & dz34m2 * ( 64.D0 )                                               
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*Qq2*gaq2*gvf2*      
     & dz34m2 * ( 64.D0 )                                               
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*Qq2*gaf2*gvq2*      
     & dz34m2 * ( 64.D0 )                                               
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*Qq2*gaf2*gaq2*      
     & dz34m2 * ( 64.D0 )                                               
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp1 )               
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp1 )               
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp1 )               
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp1 )                         
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp1 )               
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp1 )            
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp1 )            
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp1 )            
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp1 )                      
      elmat2 = elmat2 + e2*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp1 )            
      elmat2 = elmat2 + e2*p2p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*    
     & kp4m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0*kp1 )                 
      elmat2 = elmat2 + e2*p2p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*    
     & kp4m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0*kp1 )       
      elmat2 = elmat2 + e2*p2p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*    
     & kp3m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0*kp1 )                 
      elmat2 = elmat2 + e2*p2p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*    
     & kp3m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0*kp1 )       
      elmat2 = elmat2 + e2*p1p4*mf2*Qf4*Qq2*e4*dg122 * (  - 32.D0*kp2*  
     &    kp3m2 )                                                       
      elmat2 = elmat2 + e2*p1p4*mq2*Qf2*Qq4*e4*dg342 * ( 32.D0*kp3*     
     &    kp2m2 )                                                       
      elmat2 = elmat2 + e2*p1p4*Qf*Qq*Qf2*Qq2*e4*dg12*dg34 * ( 64.D0 )  
      elmat2 = elmat2 + e2*p1p4*kp4m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*dg34    
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p1p4*kp3m1*Qf4*Qq2*e4*dg122 * ( 32.D0*kp2 )  
      elmat2 = elmat2 + e2*p1p4*kp3m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*dg34    
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p1p4*kp2m1*Qf2*Qq4*e4*dg342 * ( 32.D0*kp3 )  
      elmat2 = elmat2 + e2*p1p4*kp2m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*dg34    
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p1p4*kp2m1*kp3m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34  
     &  * (  - 64.D0*p2p32 )                                            
      elmat2 = elmat2 + e2*p1p4*kp1m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*dg34    
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p1p4*kp1m1*kp4m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*  
     & dg34 * ( 32.D0*kp2 )                                             
      elmat2 = elmat2 + e2*p1p4*kp1m1*kp4m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*  
     & dg34 * (  - 32.D0*kp3 )                                          
      elmat2 = elmat2 + e2*p1p4*kp1m1*kp4m1*Qf*Qq*mq2*mf2*Qf2*Qq2*e4*   
     & dg12*dg34 * (  - 128.D0 )                                        
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*mf2*Qf2*gvf2*gvq2*dz12m2 * (   
     &     - 32.D0*kp2*kp3m2 )                                          
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*mf2*Qf2*gaq2*gvf2*dz12m2 * (   
     &     - 32.D0*kp2*kp3m2 )                                          
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*mf2*Qf2*gaf2*gvq2*dz12m2 * (   
     &     - 32.D0*kp2*kp3m2 )                                          
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*mf2*Qf2*gaf2*gaq2*dz12m2 * (   
     &     - 32.D0*kp2*kp3m2 )                                          
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*mq2*Qq2*gvf2*gvq2*dz34m2 * (   
     &    32.D0*kp3*kp2m2 )                                             
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*mq2*Qq2*gaq2*gvf2*dz34m2 * (   
     &    32.D0*kp3*kp2m2 )                                             
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*mq2*Qq2*gaf2*gvq2*dz34m2 * (   
     &    32.D0*kp3*kp2m2 )                                             
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*mq2*Qq2*gaf2*gaq2*dz34m2 * (   
     &    32.D0*kp3*kp2m2 )                                             
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*Qf*Qq*gvf2*gvq2*dz12m2*dz34m2* 
     & gzmz2 * ( 64.D0 )                                                
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*Qf*Qq*gvf2*gvq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*Qf*Qq*gaq2*gvf2*dz12m2*dz34m2* 
     & gzmz2 * ( 64.D0 )                                                
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*Qf*Qq*gaq2*gvf2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*Qf*Qq*gaf2*gvq2*dz12m2*dz34m2* 
     & gzmz2 * ( 64.D0 )                                                
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*Qf*Qq*gaf2*gvq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*Qf*Qq*gaf2*gaq2*dz12m2*dz34m2* 
     & gzmz2 * ( 64.D0 )                                                
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*Qf*Qq*gaf2*gaq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * ( 64.D0 )                                      
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf2*gvf2*gvq2*dz12m2     
     &  * ( 32.D0*kp2 )                                                 
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf2*gaq2*gvf2*dz12m2     
     &  * ( 32.D0*kp2 )                                                 
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf2*gaf2*gvq2*dz12m2     
     &  * ( 32.D0*kp2 )                                                 
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf2*gaf2*gaq2*dz12m2     
     &  * ( 32.D0*kp2 )                                                 
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qq2*gvf2*gvq2*dz34m2     
     &  * ( 32.D0*kp3 )                                                 
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qq2*gaq2*gvf2*dz34m2     
     &  * ( 32.D0*kp3 )                                                 
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qq2*gaf2*gvq2*dz34m2     
     &  * ( 32.D0*kp3 )                                                 
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qq2*gaf2*gaq2*dz34m2     
     &  * ( 32.D0*kp3 )                                                 
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gvf2*gvq2*   
     & dz12m2*dz34m2*gzmz2 * (  - 64.D0*p2p32 )                         
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gvf2*gvq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0*p2p32 )               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaq2*gvf2*   
     & dz12m2*dz34m2*gzmz2 * (  - 64.D0*p2p32 )                         
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaq2*gvf2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0*p2p32 )               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*gvq2*   
     & dz12m2*dz34m2*gzmz2 * (  - 64.D0*p2p32 )                         
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*gvq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0*p2p32 )               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*gaq2*   
     & dz12m2*dz34m2*gzmz2 * (  - 64.D0*p2p32 )                         
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*gaq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0*p2p32 )               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*gvf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp2 )                         
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*gvf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 )               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*gaq2*    
     & gvf2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp2 )                         
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*gaq2*    
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 )               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*gaf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp2 )                      
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*gaf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 )            
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*gaf2*    
     & gaq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp2 )                      
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*gaf2*    
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 )            
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*gvf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*gvf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp3 )            
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*gaq2*    
     & gvf2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*gaq2*    
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp3 )               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*gaf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*gaf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp3 )            
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*gaf2*    
     & gaq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*gaf2*    
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp3 )               
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*mf2*gvf2 
     & *gvq2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                        
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*mf2*gvf2 
     & *gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )              
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*mf2*gaq2 
     & *gvf2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                           
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*mf2*gaq2 
     & *gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                 
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*mf2*gaf2 
     & *gvq2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                           
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*mf2*gaf2 
     & *gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                 
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*mf2*gaf2 
     & *gaq2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                        
      elmat2 = elmat2 + e2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*mf2*gaf2 
     & *gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )              
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*mf2*Qf2*dz12m2 
     &  * (  - 128.D0*kp2*kp3m2 )                                       
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*mq2*Qq2*dz34m2 
     &  * ( 128.D0*kp3*kp2m2 )                                          
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*Qf*Qq*dz12m2*  
     & dz34m2*gzmz2 * ( 256.D0 )                                        
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*Qf*Qq*dz12m2*  
     & dz34m2*s12mmz2*s34mmz2 * ( 256.D0 )                              
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf2*     
     & dz12m2 * ( 128.D0*kp2 )                                          
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf*Qq*   
     & mf2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                          
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf*Qq*   
     & mf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qq2*     
     & dz34m2 * ( 128.D0*kp3 )                                          
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf*Qq*   
     & mq2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                             
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf*Qq*   
     & mq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                   
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*kp3m1*Qf 
     & *Qq*dz12m2*dz34m2*gzmz2 * (  - 256.D0*p2p32 )                    
      elmat2 = elmat2 + e2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*kp3m1*Qf 
     & *Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 256.D0*p2p32 )          
      elmat2 = elmat2 + e2*p1p4*p3p4*kp4m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p1p4*p3p4*kp3m1*kp4m1*Qf4*Qq2*e4*dg122 * (   
     &    32.D0*kp2 )                                                   
      elmat2 = elmat2 + e2*p1p4*p3p4*kp1m1*kp4m1*Qf*Qq*mq2*Qf2*Qq2*e4*  
     & dg12*dg34 * (  - 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gvf2*gvq2 
     & *dz12m2 * ( 32.D0*kp2 )                                          
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaq2*gvf2 
     & *dz12m2 * ( 32.D0*kp2 )                                          
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2*gvq2 
     & *dz12m2 * ( 32.D0*kp2 )                                          
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2*gaq2 
     & *dz12m2 * ( 32.D0*kp2 )                                          
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*    
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*    
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*    
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*    
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*    
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*    
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*    
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p3p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mq2*    
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p4*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                              
      elmat2 = elmat2 + e2*p1p4*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                    
      elmat2 = elmat2 + e2*p1p4*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*    
     & kp4m1*Qf2*dz12m2 * ( 128.D0*kp2 )                                
      elmat2 = elmat2 + e2*p1p4*p2p4*kp4m1*Qf4*Qq2*e4*dg122 * (  - 64.D0
     &     )                                                            
      elmat2 = elmat2 + e2*p1p4*p2p4*kp2m1*kp4m1*Qf*Qq*Qf2*Qq2*e4*dg12* 
     & dg34 * (  - 32.D0*kp3 )                                          
      elmat2 = elmat2 + e2*p1p4*p2p4*kp1m1*kp4m1*Qf*Qq*Qf2*Qq2*e4*dg12* 
     & dg34 * ( 32.D0*kp3 )                                             
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp4m1*Qf2*gvf2*gvq2*      
     & dz12m2 * (  - 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp4m1*Qf2*gaq2*gvf2*      
     & dz12m2 * (  - 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp4m1*Qf2*gaf2*gvq2*      
     & dz12m2 * (  - 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp4m1*Qf2*gaf2*gaq2*      
     & dz12m2 * (  - 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp3 )            
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp3 )            
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp3 )            
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp3 )                      
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp3 )            
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp3 )               
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp3 )               
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp3 )               
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp3 )                         
      elmat2 = elmat2 + e2*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp3 )               
      elmat2 = elmat2 + e2*p1p4*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*    
     & kp4m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0*kp3 )              
      elmat2 = elmat2 + e2*p1p4*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*    
     & kp4m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0*kp3 )    
      elmat2 = elmat2 + e2*p1p4*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp4m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0*kp3 )              
      elmat2 = elmat2 + e2*p1p4*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp4m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0*kp3 )    
      elmat2 = elmat2 + e2*p1p4*p2p3*mf2*Qf4*Qq2*e4*dg122 * (  - 32.D0* 
     &    kp4m2 - 32.D0*kp3m2 )                                         
      elmat2 = elmat2 + e2*p1p4*p2p3*mq2*Qf2*Qq4*e4*dg342 * (  - 32.D0* 
     &    kp2m2 - 32.D0*kp1m2 )                                         
      elmat2 = elmat2 + e2*p1p4*p2p3*kp4m1*Qf4*Qq2*e4*dg122 * ( 32.D0 ) 
      elmat2 = elmat2 + e2*p1p4*p2p3*kp4m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p1p4*p2p3*kp3m1*Qf4*Qq2*e4*dg122 * ( 32.D0 ) 
      elmat2 = elmat2 + e2*p1p4*p2p3*kp3m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p1p4*p2p3*kp2m1*Qf2*Qq4*e4*dg342 * (  - 32.D0
     &     )                                                            
      elmat2 = elmat2 + e2*p1p4*p2p3*kp2m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p1p4*p2p3*kp1m1*Qf2*Qq4*e4*dg342 * (  - 32.D0
     &     )                                                            
      elmat2 = elmat2 + e2*p1p4*p2p3*kp1m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*mf2*Qf2*gvf2*gvq2*dz12m2  
     &  * (  - 32.D0*kp4m2 - 32.D0*kp3m2 )                              
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*mf2*Qf2*gaq2*gvf2*dz12m2  
     &  * (  - 32.D0*kp4m2 - 32.D0*kp3m2 )                              
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*mf2*Qf2*gaf2*gvq2*dz12m2  
     &  * (  - 32.D0*kp4m2 - 32.D0*kp3m2 )                              
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*mf2*Qf2*gaf2*gaq2*dz12m2  
     &  * (  - 32.D0*kp4m2 - 32.D0*kp3m2 )                              
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*mq2*Qq2*gvf2*gvq2*dz34m2  
     &  * (  - 32.D0*kp2m2 - 32.D0*kp1m2 )                              
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*mq2*Qq2*gaq2*gvf2*dz34m2  
     &  * (  - 32.D0*kp2m2 - 32.D0*kp1m2 )                              
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*mq2*Qq2*gaf2*gvq2*dz34m2  
     &  * (  - 32.D0*kp2m2 - 32.D0*kp1m2 )                              
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*mq2*Qq2*gaf2*gaq2*dz34m2  
     &  * (  - 32.D0*kp2m2 - 32.D0*kp1m2 )                              
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf2*gvf2*gvq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf2*gaq2*gvf2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf2*gaf2*gvq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf2*gaf2*gaq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp4m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf2*gvf2*gvq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf2*gaq2*gvf2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf2*gaf2*gvq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf2*gaf2*gaq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp3m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qq2*gvf2*gvq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qq2*gaq2*gvf2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qq2*gaf2*gvq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qq2*gaf2*gaq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qq2*gvf2*gvq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qq2*gaq2*gvf2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qq2*gaf2*gvq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qq2*gaf2*gaq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p4*p2p3*gs2sctw4*kp1m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*mf2*Qf2*  
     & dz12m2 * (  - 128.D0*kp4m2 - 128.D0*kp3m2 )                      
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*mq2*Qq2*  
     & dz34m2 * (  - 128.D0*kp2m2 - 128.D0*kp1m2 )                      
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf2 
     & *dz12m2 * ( 128.D0 )                                             
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                           
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                 
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf2 
     & *dz12m2 * ( 128.D0 )                                             
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                           
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                 
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qq2 
     & *dz34m2 * (  - 128.D0 )                                          
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                              
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                    
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qq2 
     & *dz34m2 * (  - 128.D0 )                                          
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                              
      elmat2 = elmat2 + e2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                    
      elmat2 = elmat2 + e2*p1p4*p2p3*p3p4*kp3m1*kp4m1*Qf4*Qq2*e4*dg122  
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p1p4*p2p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gvf2 
     & *gvq2*dz12m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p4*p2p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaq2 
     & *gvf2*dz12m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p4*p2p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2 
     & *gvq2*dz12m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p4*p2p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2 
     & *gaq2*dz12m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p4*p2p3*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp3m1*kp4m1*Qf2*dz12m2 * ( 256.D0 )                              
      elmat2 = elmat2 + e2*p1p4*p2p3*p2p4*kp2m1*kp4m1*Qf*Qq*Qf2*Qq2*e4* 
     & dg12*dg34 * ( 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p4*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*   
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*   
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p4*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*   
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*   
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p4*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*   
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*   
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p4*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*   
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p4*p2p3*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*   
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p4*p2p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp2m1*kp4m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * ( 256.D0 )               
      elmat2 = elmat2 + e2*p1p4*p2p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp2m1*kp4m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 256.D0 )     
      elmat2 = elmat2 + e2*p1p3*mf2*Qf4*Qq2*e4*dg122 * (  - 32.D0*kp2*  
     &    kp4m2 )                                                       
      elmat2 = elmat2 + e2*p1p3*mq2*Qf2*Qq4*e4*dg342 * ( 32.D0*kp4*     
     &    kp2m2 )                                                       
      elmat2 = elmat2 + e2*p1p3*Qf*Qq*Qf2*Qq2*e4*dg12*dg34 * (  - 64.D0 
     &     )                                                            
      elmat2 = elmat2 + e2*p1p3*kp4m1*Qf4*Qq2*e4*dg122 * ( 32.D0*kp2 )  
      elmat2 = elmat2 + e2*p1p3*kp4m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*dg34    
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p1p3*kp3m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*dg34    
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p1p3*kp2m1*Qf2*Qq4*e4*dg342 * ( 32.D0*kp4 )  
      elmat2 = elmat2 + e2*p1p3*kp2m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*dg34    
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p1p3*kp2m1*kp4m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34  
     &  * ( 64.D0*p2p42 )                                               
      elmat2 = elmat2 + e2*p1p3*kp1m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*dg34    
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p1p3*kp1m1*kp3m1*Qf*Qq*mf2*Qf2*Qq2*e4*dg12*  
     & dg34 * (  - 32.D0*kp2 )                                          
      elmat2 = elmat2 + e2*p1p3*kp1m1*kp3m1*Qf*Qq*mq2*Qf2*Qq2*e4*dg12*  
     & dg34 * ( 32.D0*kp4 )                                             
      elmat2 = elmat2 + e2*p1p3*kp1m1*kp3m1*Qf*Qq*mq2*mf2*Qf2*Qq2*e4*   
     & dg12*dg34 * ( 128.D0 )                                           
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*mf2*Qf2*gvf2*gvq2*dz12m2 * (   
     &     - 32.D0*kp2*kp4m2 )                                          
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*mf2*Qf2*gaq2*gvf2*dz12m2 * (   
     &     - 32.D0*kp2*kp4m2 )                                          
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*mf2*Qf2*gaf2*gvq2*dz12m2 * (   
     &     - 32.D0*kp2*kp4m2 )                                          
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*mf2*Qf2*gaf2*gaq2*dz12m2 * (   
     &     - 32.D0*kp2*kp4m2 )                                          
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*mq2*Qq2*gvf2*gvq2*dz34m2 * (   
     &    32.D0*kp4*kp2m2 )                                             
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*mq2*Qq2*gaq2*gvf2*dz34m2 * (   
     &    32.D0*kp4*kp2m2 )                                             
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*mq2*Qq2*gaf2*gvq2*dz34m2 * (   
     &    32.D0*kp4*kp2m2 )                                             
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*mq2*Qq2*gaf2*gaq2*dz34m2 * (   
     &    32.D0*kp4*kp2m2 )                                             
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*Qf*Qq*gvf2*gvq2*dz12m2*dz34m2* 
     & gzmz2 * (  - 64.D0 )                                             
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*Qf*Qq*gvf2*gvq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*Qf*Qq*gaq2*gvf2*dz12m2*dz34m2* 
     & gzmz2 * (  - 64.D0 )                                             
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*Qf*Qq*gaq2*gvf2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*Qf*Qq*gaf2*gvq2*dz12m2*dz34m2* 
     & gzmz2 * (  - 64.D0 )                                             
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*Qf*Qq*gaf2*gvq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*Qf*Qq*gaf2*gaq2*dz12m2*dz34m2* 
     & gzmz2 * (  - 64.D0 )                                             
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*Qf*Qq*gaf2*gaq2*dz12m2*dz34m2* 
     & s12mmz2*s34mmz2 * (  - 64.D0 )                                   
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf2*gvf2*gvq2*dz12m2     
     &  * ( 32.D0*kp2 )                                                 
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf2*gaq2*gvf2*dz12m2     
     &  * ( 32.D0*kp2 )                                                 
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf2*gaf2*gvq2*dz12m2     
     &  * ( 32.D0*kp2 )                                                 
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf2*gaf2*gaq2*dz12m2     
     &  * ( 32.D0*kp2 )                                                 
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp4m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp3m1*Qf*Qq*mf2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qq2*gvf2*gvq2*dz34m2     
     &  * ( 32.D0*kp4 )                                                 
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qq2*gaq2*gvf2*dz34m2     
     &  * ( 32.D0*kp4 )                                                 
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qq2*gaf2*gvq2*dz34m2     
     &  * ( 32.D0*kp4 )                                                 
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qq2*gaf2*gaq2*dz34m2     
     &  * ( 32.D0*kp4 )                                                 
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gvf2*gvq2*   
     & dz12m2*dz34m2*gzmz2 * ( 64.D0*p2p42 )                            
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gvf2*gvq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0*p2p42 )                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaq2*gvf2*   
     & dz12m2*dz34m2*gzmz2 * ( 64.D0*p2p42 )                            
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaq2*gvf2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0*p2p42 )                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*gvq2*   
     & dz12m2*dz34m2*gzmz2 * ( 64.D0*p2p42 )                            
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*gvq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0*p2p42 )                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*gaq2*   
     & dz12m2*dz34m2*gzmz2 * ( 64.D0*p2p42 )                            
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp2m1*kp4m1*Qf*Qq*gaf2*gaq2*   
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0*p2p42 )                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gvf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaq2*gvf2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gvq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*Qf*Qq*mq2*gaf2*gaq2*     
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*gvf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp2 )                      
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*gvf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 )            
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*gaq2*    
     & gvf2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp2 )                      
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*gaq2*    
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 )            
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*gaf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp2 )                         
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*gaf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 )               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*gaf2*    
     & gaq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp2 )                         
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*gaf2*    
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 )               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*gvf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp4 )                         
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*gvf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 )               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*gaq2*    
     & gvf2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp4 )                      
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*gaq2*    
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 )            
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*gaf2*    
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp4 )                         
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*gaf2*    
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 )               
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*gaf2*    
     & gaq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp4 )                      
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*gaf2*    
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 )            
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*mf2*gvf2 
     & *gvq2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                           
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*mf2*gvf2 
     & *gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                 
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*mf2*gaq2 
     & *gvf2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                        
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*mf2*gaq2 
     & *gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )              
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*mf2*gaf2 
     & *gvq2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                        
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*mf2*gaf2 
     & *gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )              
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*mf2*gaf2 
     & *gaq2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                           
      elmat2 = elmat2 + e2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*mf2*gaf2 
     & *gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                 
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*mf2*Qf2*dz12m2 
     &  * ( 128.D0*kp2*kp4m2 )                                          
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*mq2*Qq2*dz34m2 
     &  * (  - 128.D0*kp4*kp2m2 )                                       
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*Qf*Qq*dz12m2*  
     & dz34m2*gzmz2 * ( 256.D0 )                                        
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*Qf*Qq*dz12m2*  
     & dz34m2*s12mmz2*s34mmz2 * ( 256.D0 )                              
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf2*     
     & dz12m2 * (  - 128.D0*kp2 )                                       
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf*Qq*   
     & mf2*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                          
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf*Qq*   
     & mf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qq2*     
     & dz34m2 * (  - 128.D0*kp4 )                                       
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf*Qq*   
     & mq2*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                             
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf*Qq*   
     & mq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                   
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*kp4m1*Qf 
     & *Qq*dz12m2*dz34m2*gzmz2 * (  - 256.D0*p2p42 )                    
      elmat2 = elmat2 + e2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*kp4m1*Qf 
     & *Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 256.D0*p2p42 )          
      elmat2 = elmat2 + e2*p1p3*p3p4*kp3m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p1p3*p3p4*kp3m1*kp4m1*Qf4*Qq2*e4*dg122 * (   
     &    32.D0*kp2 )                                                   
      elmat2 = elmat2 + e2*p1p3*p3p4*kp1m1*kp3m1*Qf*Qq*mq2*Qf2*Qq2*e4*  
     & dg12*dg34 * ( 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gvf2*gvq2 
     & *dz12m2 * ( 32.D0*kp2 )                                          
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaq2*gvf2 
     & *dz12m2 * ( 32.D0*kp2 )                                          
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2*gvq2 
     & *dz12m2 * ( 32.D0*kp2 )                                          
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2*gaq2 
     & *dz12m2 * ( 32.D0*kp2 )                                          
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*    
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*    
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*    
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*    
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*    
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*    
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*    
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p3p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mq2*    
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p3*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                              
      elmat2 = elmat2 + e2*p1p3*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                    
      elmat2 = elmat2 + e2*p1p3*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*    
     & kp4m1*Qf2*dz12m2 * (  - 128.D0*kp2 )                             
      elmat2 = elmat2 + e2*p1p3*p2p4*mf2*Qf4*Qq2*e4*dg122 * (  - 32.D0* 
     &    kp4m2 - 32.D0*kp3m2 )                                         
      elmat2 = elmat2 + e2*p1p3*p2p4*mq2*Qf2*Qq4*e4*dg342 * (  - 32.D0* 
     &    kp2m2 - 32.D0*kp1m2 )                                         
      elmat2 = elmat2 + e2*p1p3*p2p4*kp4m1*Qf4*Qq2*e4*dg122 * ( 32.D0 ) 
      elmat2 = elmat2 + e2*p1p3*p2p4*kp4m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p1p3*p2p4*kp3m1*Qf4*Qq2*e4*dg122 * ( 32.D0 ) 
      elmat2 = elmat2 + e2*p1p3*p2p4*kp3m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p1p3*p2p4*kp2m1*Qf2*Qq4*e4*dg342 * (  - 32.D0
     &     )                                                            
      elmat2 = elmat2 + e2*p1p3*p2p4*kp2m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p1p3*p2p4*kp1m1*Qf2*Qq4*e4*dg342 * (  - 32.D0
     &     )                                                            
      elmat2 = elmat2 + e2*p1p3*p2p4*kp1m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*mf2*Qf2*gvf2*gvq2*dz12m2  
     &  * (  - 32.D0*kp4m2 - 32.D0*kp3m2 )                              
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*mf2*Qf2*gaq2*gvf2*dz12m2  
     &  * (  - 32.D0*kp4m2 - 32.D0*kp3m2 )                              
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*mf2*Qf2*gaf2*gvq2*dz12m2  
     &  * (  - 32.D0*kp4m2 - 32.D0*kp3m2 )                              
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*mf2*Qf2*gaf2*gaq2*dz12m2  
     &  * (  - 32.D0*kp4m2 - 32.D0*kp3m2 )                              
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*mq2*Qq2*gvf2*gvq2*dz34m2  
     &  * (  - 32.D0*kp2m2 - 32.D0*kp1m2 )                              
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*mq2*Qq2*gaq2*gvf2*dz34m2  
     &  * (  - 32.D0*kp2m2 - 32.D0*kp1m2 )                              
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*mq2*Qq2*gaf2*gvq2*dz34m2  
     &  * (  - 32.D0*kp2m2 - 32.D0*kp1m2 )                              
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*mq2*Qq2*gaf2*gaq2*dz34m2  
     &  * (  - 32.D0*kp2m2 - 32.D0*kp1m2 )                              
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf2*gvf2*gvq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf2*gaq2*gvf2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf2*gaf2*gvq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf2*gaf2*gaq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp4m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf2*gvf2*gvq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf2*gaq2*gvf2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf2*gaf2*gvq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf2*gaf2*gaq2*      
     & dz12m2 * ( 32.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp3m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qq2*gvf2*gvq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qq2*gaq2*gvf2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qq2*gaf2*gvq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qq2*gaf2*gaq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qq2*gvf2*gvq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qq2*gaq2*gvf2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qq2*gaf2*gvq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qq2*gaf2*gaq2*      
     & dz34m2 * (  - 32.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p3*p2p4*gs2sctw4*kp1m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*mf2*Qf2*  
     & dz12m2 * ( 128.D0*kp4m2 + 128.D0*kp3m2 )                         
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*mq2*Qq2*  
     & dz34m2 * ( 128.D0*kp2m2 + 128.D0*kp1m2 )                         
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf2 
     & *dz12m2 * (  - 128.D0 )                                          
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                           
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp4m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                 
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf2 
     & *dz12m2 * (  - 128.D0 )                                          
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                           
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp3m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                 
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qq2 
     & *dz34m2 * ( 128.D0 )                                             
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                              
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                    
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qq2 
     & *dz34m2 * ( 128.D0 )                                             
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0 )                              
      elmat2 = elmat2 + e2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0 )                    
      elmat2 = elmat2 + e2*p1p3*p2p4*p3p4*kp3m1*kp4m1*Qf4*Qq2*e4*dg122  
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p1p3*p2p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gvf2 
     & *gvq2*dz12m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p3*p2p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaq2 
     & *gvf2*dz12m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p3*p2p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2 
     & *gvq2*dz12m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p3*p2p4*p3p4*gs2sctw4*kp3m1*kp4m1*Qf2*gaf2 
     & *gaq2*dz12m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p3*p2p4*p3p4*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp3m1*kp4m1*Qf2*dz12m2 * (  - 256.D0 )                           
      elmat2 = elmat2 + e2*p1p3*p2p3*kp3m1*Qf4*Qq2*e4*dg122 * (  - 64.D0
     &     )                                                            
      elmat2 = elmat2 + e2*p1p3*p2p3*kp2m1*kp3m1*Qf*Qq*Qf2*Qq2*e4*dg12* 
     & dg34 * ( 32.D0*kp4 )                                             
      elmat2 = elmat2 + e2*p1p3*p2p3*kp1m1*kp3m1*Qf*Qq*Qf2*Qq2*e4*dg12* 
     & dg34 * (  - 32.D0*kp4 )                                          
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp3m1*Qf2*gvf2*gvq2*      
     & dz12m2 * (  - 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp3m1*Qf2*gaq2*gvf2*      
     & dz12m2 * (  - 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp3m1*Qf2*gaf2*gvq2*      
     & dz12m2 * (  - 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp3m1*Qf2*gaf2*gaq2*      
     & dz12m2 * (  - 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp4 )                         
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 )               
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp4 )                         
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 )               
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp4 )                         
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 )               
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp4 )                         
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp4 )               
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp4 )                      
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 )            
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp4 )                      
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 )            
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp4 )                      
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 )            
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp4 )                      
      elmat2 = elmat2 + e2*p1p3*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp4 )            
      elmat2 = elmat2 + e2*p1p3*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*    
     & kp3m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0*kp4 )              
      elmat2 = elmat2 + e2*p1p3*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*    
     & kp3m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0*kp4 )    
      elmat2 = elmat2 + e2*p1p3*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp3m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0*kp4 )              
      elmat2 = elmat2 + e2*p1p3*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp3m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0*kp4 )    
      elmat2 = elmat2 + e2*p1p3*p2p3*p2p4*kp2m1*kp3m1*Qf*Qq*Qf2*Qq2*e4* 
     & dg12*dg34 * (  - 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p3*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*   
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*   
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p3*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*   
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*   
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p3*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*   
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*   
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p3*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*   
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p2p3*p2p4*gs2sctw4*kp2m1*kp3m1*Qf*Qq*   
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p3*p2p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp2m1*kp3m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * ( 256.D0 )               
      elmat2 = elmat2 + e2*p1p3*p2p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp2m1*kp3m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 256.D0 )     
      elmat2 = elmat2 + e2*p1p3*p1p4*kp1m1*Qf2*Qq4*e4*dg342 * ( 64.D0 ) 
      elmat2 = elmat2 + e2*p1p3*p1p4*kp1m1*kp4m1*Qf*Qq*Qf2*Qq2*e4*dg12* 
     & dg34 * (  - 32.D0*kp2 )                                          
      elmat2 = elmat2 + e2*p1p3*p1p4*kp1m1*kp3m1*Qf*Qq*Qf2*Qq2*e4*dg12* 
     & dg34 * ( 32.D0*kp2 )                                             
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*Qq2*gvf2*gvq2*      
     & dz34m2 * ( 64.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*Qq2*gaq2*gvf2*      
     & dz34m2 * ( 64.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*Qq2*gaf2*gvq2*      
     & dz34m2 * ( 64.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*Qq2*gaf2*gaq2*      
     & dz34m2 * ( 64.D0 )                                               
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp2 )                      
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 )            
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp2 )                      
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 )            
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp2 )                      
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 )            
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*gzmz2 * (  - 32.D0*kp2 )                      
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0*kp2 )            
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp2 )                         
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gvf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 )               
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp2 )                         
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaq2*   
     & gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 )               
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp2 )                         
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*   
     & gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 )               
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*gzmz2 * ( 32.D0*kp2 )                         
      elmat2 = elmat2 + e2*p1p3*p1p4*gs2sctw4*kp1m1*kp3m1*Qf*Qq*gaf2*   
     & gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0*kp2 )               
      elmat2 = elmat2 + e2*p1p3*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp4m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0*kp2 )                 
      elmat2 = elmat2 + e2*p1p3*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp4m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0*kp2 )       
      elmat2 = elmat2 + e2*p1p3*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp3m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * ( 128.D0*kp2 )                 
      elmat2 = elmat2 + e2*p1p3*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp3m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 128.D0*kp2 )       
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p4*kp1m1*kp4m1*Qf*Qq*Qf2*Qq2*e4* 
     & dg12*dg34 * (  - 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*   
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*   
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*   
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*   
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*   
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*   
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*   
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*   
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp1m1*kp4m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * ( 256.D0 )               
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp1m1*kp4m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 256.D0 )     
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p3*kp1m1*kp3m1*Qf*Qq*Qf2*Qq2*e4* 
     & dg12*dg34 * ( 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*   
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*   
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*   
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*   
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*   
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*   
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*   
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*   
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp1m1*kp3m1*Qf*Qq*dz12m2*dz34m2*gzmz2 * ( 256.D0 )               
      elmat2 = elmat2 + e2*p1p3*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp1m1*kp3m1*Qf*Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 256.D0 )     
      elmat2 = elmat2 + e2*p1p2*mf4*Qf4*Qq2*e4*dg122 * (  - 32.D0*kp4m2 
     &     - 32.D0*kp3m2 )                                              
      elmat2 = elmat2 + e2*p1p2*mq2*mf2*Qf2*Qq4*e4*dg342 * (  - 32.D0*  
     &    kp2m2 - 32.D0*kp1m2 )                                         
      elmat2 = elmat2 + e2*p1p2*kp2m1*mf2*Qf2*Qq4*e4*dg342 * (  - 64.D0 
     &     )                                                            
      elmat2 = elmat2 + e2*p1p2*kp1m1*mf2*Qf2*Qq4*e4*dg342 * (  - 64.D0 
     &     )                                                            
      elmat2 = elmat2 + e2*p1p2*kp1m1*kp2m1*mq2*mf2*Qf2*Qq4*e4*dg342    
     &  * ( 128.D0 )                                                    
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*mf4*Qf2*gvf2*gvq2*dz12m2 * (   
     &     - 32.D0*kp4m2 - 32.D0*kp3m2 )                                
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*mf4*Qf2*gaq2*gvf2*dz12m2 * (   
     &     - 32.D0*kp4m2 - 32.D0*kp3m2 )                                
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*mf4*Qf2*gaf2*gvq2*dz12m2 * (   
     &    32.D0*kp4m2 + 32.D0*kp3m2 )                                   
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*mf4*Qf2*gaf2*gaq2*dz12m2 * (   
     &    32.D0*kp4m2 + 32.D0*kp3m2 )                                   
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*mq2*mf2*Qq2*gvf2*gvq2*dz34m2   
     &  * (  - 32.D0*kp2m2 - 32.D0*kp1m2 )                              
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*mq2*mf2*Qq2*gaq2*gvf2*dz34m2   
     &  * (  - 32.D0*kp2m2 - 32.D0*kp1m2 )                              
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*mq2*mf2*Qq2*gaf2*gvq2*dz34m2   
     &  * ( 32.D0*kp2m2 + 32.D0*kp1m2 )                                 
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*mq2*mf2*Qq2*gaf2*gaq2*dz34m2   
     &  * ( 32.D0*kp2m2 + 32.D0*kp1m2 )                                 
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp2m1*mf2*Qq2*gvf2*gvq2*dz34m2 
     &  * (  - 64.D0 )                                                  
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp2m1*mf2*Qq2*gaq2*gvf2*dz34m2 
     &  * (  - 64.D0 )                                                  
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp2m1*mf2*Qq2*gaf2*gvq2*dz34m2 
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp2m1*mf2*Qq2*gaf2*gaq2*dz34m2 
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp1m1*mf2*Qq2*gvf2*gvq2*dz34m2 
     &  * (  - 64.D0 )                                                  
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp1m1*mf2*Qq2*gaq2*gvf2*dz34m2 
     &  * (  - 64.D0 )                                                  
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp1m1*mf2*Qq2*gaf2*gvq2*dz34m2 
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp1m1*mf2*Qq2*gaf2*gaq2*dz34m2 
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp1m1*kp2m1*mq2*mf2*Qq2*gvf2*  
     & gvq2*dz34m2 * ( 128.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp1m1*kp2m1*mq2*mf2*Qq2*gaq2*  
     & gvf2*dz34m2 * (  - 128.D0 )                                      
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp1m1*kp2m1*mq2*mf2*Qq2*gaf2*  
     & gvq2*dz34m2 * (  - 128.D0 )                                      
      elmat2 = elmat2 + e2*p1p2*gs2sctw4*kp1m1*kp2m1*mq2*mf2*Qq2*gaf2*  
     & gaq2*dz34m2 * ( 128.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*p3p4*kp3m1*kp4m1*mf2*Qf4*Qq2*e4*dg122   
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p1p2*p3p4*kp1m1*kp2m1*mq2*Qf2*Qq4*e4*dg342   
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p1p2*p3p4*gs2sctw4*kp3m1*kp4m1*mf2*Qf2*gvf2* 
     & gvq2*dz12m2 * ( 64.D0 )                                          
      elmat2 = elmat2 + e2*p1p2*p3p4*gs2sctw4*kp3m1*kp4m1*mf2*Qf2*gaq2* 
     & gvf2*dz12m2 * ( 64.D0 )                                          
      elmat2 = elmat2 + e2*p1p2*p3p4*gs2sctw4*kp3m1*kp4m1*mf2*Qf2*gaf2* 
     & gvq2*dz12m2 * (  - 64.D0 )                                       
      elmat2 = elmat2 + e2*p1p2*p3p4*gs2sctw4*kp3m1*kp4m1*mf2*Qf2*gaf2* 
     & gaq2*dz12m2 * (  - 64.D0 )                                       
      elmat2 = elmat2 + e2*p1p2*p3p4*gs2sctw4*kp1m1*kp2m1*mq2*Qq2*gvf2* 
     & gvq2*dz34m2 * ( 64.D0 )                                          
      elmat2 = elmat2 + e2*p1p2*p3p4*gs2sctw4*kp1m1*kp2m1*mq2*Qq2*gaq2* 
     & gvf2*dz34m2 * (  - 64.D0 )                                       
      elmat2 = elmat2 + e2*p1p2*p3p4*gs2sctw4*kp1m1*kp2m1*mq2*Qq2*gaf2* 
     & gvq2*dz34m2 * ( 64.D0 )                                          
      elmat2 = elmat2 + e2*p1p2*p3p4*gs2sctw4*kp1m1*kp2m1*mq2*Qq2*gaf2* 
     & gaq2*dz34m2 * (  - 64.D0 )                                       
      elmat2 = elmat2 + e2*p1p2*p2p4*kp2m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p1p2*p2p4*kp2m1*kp4m1*Qf*Qq*mf2*Qf2*Qq2*e4*  
     & dg12*dg34 * ( 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p2*p2p4*kp1m1*kp2m1*Qf2*Qq4*e4*dg342 * (   
     &     - 32.D0*kp3 )                                                
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*    
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*    
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*    
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*    
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*    
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*    
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*    
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp2m1*kp4m1*Qf*Qq*mf2*    
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp1m1*kp2m1*Qq2*gvf2*gvq2 
     & *dz34m2 * (  - 32.D0*kp3 )                                       
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp1m1*kp2m1*Qq2*gaq2*gvf2 
     & *dz34m2 * (  - 32.D0*kp3 )                                       
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2*gvq2 
     & *dz34m2 * (  - 32.D0*kp3 )                                       
      elmat2 = elmat2 + e2*p1p2*p2p4*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2*gaq2 
     & *dz34m2 * (  - 32.D0*kp3 )                                       
      elmat2 = elmat2 + e2*p1p2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                           
      elmat2 = elmat2 + e2*p1p2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                 
      elmat2 = elmat2 + e2*p1p2*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp2m1*Qq2*dz34m2 * ( 128.D0*kp3 )                                
      elmat2 = elmat2 + e2*p1p2*p2p3*kp2m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p1p2*p2p3*kp2m1*kp3m1*Qf*Qq*mf2*Qf2*Qq2*e4*  
     & dg12*dg34 * (  - 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*p2p3*kp1m1*kp2m1*Qf2*Qq4*e4*dg342 * (   
     &     - 32.D0*kp4 )                                                
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*    
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*    
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*    
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*    
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*    
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*    
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*    
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp2m1*kp3m1*Qf*Qq*mf2*    
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp1m1*kp2m1*Qq2*gvf2*gvq2 
     & *dz34m2 * (  - 32.D0*kp4 )                                       
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp1m1*kp2m1*Qq2*gaq2*gvf2 
     & *dz34m2 * (  - 32.D0*kp4 )                                       
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2*gvq2 
     & *dz34m2 * (  - 32.D0*kp4 )                                       
      elmat2 = elmat2 + e2*p1p2*p2p3*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2*gaq2 
     & *dz34m2 * (  - 32.D0*kp4 )                                       
      elmat2 = elmat2 + e2*p1p2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                           
      elmat2 = elmat2 + e2*p1p2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp2m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                 
      elmat2 = elmat2 + e2*p1p2*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp2m1*Qq2*dz34m2 * (  - 128.D0*kp4 )                             
      elmat2 = elmat2 + e2*p1p2*p1p4*kp1m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * (  - 32.D0 )                                                  
      elmat2 = elmat2 + e2*p1p2*p1p4*kp1m1*kp4m1*Qf*Qq*mf2*Qf2*Qq2*e4*  
     & dg12*dg34 * (  - 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*p1p4*kp1m1*kp2m1*Qf2*Qq4*e4*dg342 * (   
     &     - 32.D0*kp3 )                                                
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * (  - 32.D0 )                               
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 32.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*    
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*    
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*    
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*    
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*    
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*    
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*    
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp4m1*Qf*Qq*mf2*    
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp2m1*Qq2*gvf2*gvq2 
     & *dz34m2 * (  - 32.D0*kp3 )                                       
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp2m1*Qq2*gaq2*gvf2 
     & *dz34m2 * (  - 32.D0*kp3 )                                       
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2*gvq2 
     & *dz34m2 * (  - 32.D0*kp3 )                                       
      elmat2 = elmat2 + e2*p1p2*p1p4*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2*gaq2 
     & *dz34m2 * (  - 32.D0*kp3 )                                       
      elmat2 = elmat2 + e2*p1p2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                           
      elmat2 = elmat2 + e2*p1p2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                 
      elmat2 = elmat2 + e2*p1p2*p1p4*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp2m1*Qq2*dz34m2 * (  - 128.D0*kp3 )                             
      elmat2 = elmat2 + e2*p1p2*p1p4*p2p3*kp1m1*kp2m1*Qf2*Qq4*e4*dg342  
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p1p2*p1p4*p2p3*gs2sctw4*kp1m1*kp2m1*Qq2*gvf2 
     & *gvq2*dz34m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*p1p4*p2p3*gs2sctw4*kp1m1*kp2m1*Qq2*gaq2 
     & *gvf2*dz34m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*p1p4*p2p3*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2 
     & *gvq2*dz34m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*p1p4*p2p3*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2 
     & *gaq2*dz34m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*p1p4*p2p3*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp1m1*kp2m1*Qq2*dz34m2 * ( 256.D0 )                              
      elmat2 = elmat2 + e2*p1p2*p1p3*kp1m1*Qf*Qq*Qf2*Qq2*e4*dg12*dg34   
     &  * ( 32.D0 )                                                     
      elmat2 = elmat2 + e2*p1p2*p1p3*kp1m1*kp3m1*Qf*Qq*mf2*Qf2*Qq2*e4*  
     & dg12*dg34 * ( 64.D0 )                                            
      elmat2 = elmat2 + e2*p1p2*p1p3*kp1m1*kp2m1*Qf2*Qq4*e4*dg342 * (   
     &     - 32.D0*kp4 )                                                
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*Qf*Qq*gvf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*Qf*Qq*gaq2*gvf2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*Qf*Qq*gaf2*gvq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*gzmz2 * ( 32.D0 )                                  
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*Qf*Qq*gaf2*gaq2*    
     & dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 32.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*    
     & gvf2*gvq2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*    
     & gvf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*    
     & gaq2*gvf2*dz12m2*dz34m2*gzmz2 * ( 64.D0 )                        
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*    
     & gaq2*gvf2*dz12m2*dz34m2*s12mmz2*s34mmz2 * ( 64.D0 )              
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*    
     & gaf2*gvq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*    
     & gaf2*gvq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*    
     & gaf2*gaq2*dz12m2*dz34m2*gzmz2 * (  - 64.D0 )                     
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp3m1*Qf*Qq*mf2*    
     & gaf2*gaq2*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 64.D0 )           
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp2m1*Qq2*gvf2*gvq2 
     & *dz34m2 * (  - 32.D0*kp4 )                                       
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp2m1*Qq2*gaq2*gvf2 
     & *dz34m2 * (  - 32.D0*kp4 )                                       
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2*gvq2 
     & *dz34m2 * (  - 32.D0*kp4 )                                       
      elmat2 = elmat2 + e2*p1p2*p1p3*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2*gaq2 
     & *dz34m2 * (  - 32.D0*kp4 )                                       
      elmat2 = elmat2 + e2*p1p2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf* 
     & Qq*dz12m2*dz34m2*gzmz2 * (  - 128.D0 )                           
      elmat2 = elmat2 + e2*p1p2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*Qf* 
     & Qq*dz12m2*dz34m2*s12mmz2*s34mmz2 * (  - 128.D0 )                 
      elmat2 = elmat2 + e2*p1p2*p1p3*gaf*gaq*gvf*gvq*gs2sctw4*kp1m1*    
     & kp2m1*Qq2*dz34m2 * ( 128.D0*kp4 )                                
      elmat2 = elmat2 + e2*p1p2*p1p3*p2p4*kp1m1*kp2m1*Qf2*Qq4*e4*dg342  
     &  * ( 64.D0 )                                                     
      elmat2 = elmat2 + e2*p1p2*p1p3*p2p4*gs2sctw4*kp1m1*kp2m1*Qq2*gvf2 
     & *gvq2*dz34m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*p1p3*p2p4*gs2sctw4*kp1m1*kp2m1*Qq2*gaq2 
     & *gvf2*dz34m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*p1p3*p2p4*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2 
     & *gvq2*dz34m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*p1p3*p2p4*gs2sctw4*kp1m1*kp2m1*Qq2*gaf2 
     & *gaq2*dz34m2 * ( 64.D0 )                                         
      elmat2 = elmat2 + e2*p1p2*p1p3*p2p4*gaf*gaq*gvf*gvq*gs2sctw4*     
     & kp1m1*kp2m1*Qq2*dz34m2 * (  - 256.D0 )                           
                                                                        
                                                                        
 !***** number of lines = 3076
