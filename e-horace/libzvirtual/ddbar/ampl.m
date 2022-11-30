
f[1] =         ((8*SW2)/(3) - (SW2)/(SW2) - 
		(4*SW2*SW2)/(3))/(-MZ2+gzmz + S);


f[2] =         ((2*SW2) - (4*SW2*SW2)/(3))/(-MZ2+gzmz + S);


f[3] =         ((2*SW2)/(3) - (4*SW2*SW2)/(3))/
(-MZ2+gzmz + S);


f[4] =        -(4*SW2*SW2)/(3)/(-MZ2+gzmz + S);

Do[ g[i] = Simplify[ D[ f[i], SW2]];Print["i=",i,"  ",g[i],"\n"] , {i,1,4}];
