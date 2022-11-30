f[1] = ((-10*SW2*Pi)/(3) + (SW2*Pi)/(SW2) + 
        (8*SW2*Pi*SW2)/(3))/(-MZ2+gzmz + S);

f[2] = ((-2*SW2*Pi) + (8*SW2*Pi*SW2)/(3))/(-MZ2+gzmz + S) ;

f[3] = ((-4*SW2*Pi)/(3) + (8*SW2*Pi*SW2)/(3))/(-MZ2+gzmz + S) ;

f[4] =  (8*SW2*Pi*SW2)/(3)/(-MZ2+gzmz + S) ;


Do[ g[i] = Simplify[ D[ f[i]/Pi, SW2]];Print["i=",i,"  ",g[i],"\n"] , {i,1,4}];
