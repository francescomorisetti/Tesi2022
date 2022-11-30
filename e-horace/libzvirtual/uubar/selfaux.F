
#include "prefix.h"

        double complex function Agg(SS)
	implicit character (a-s,u-z)
	implicit double complex (t)
        double precision SS


#include "vars.h"


      Agg= (Alfa*SS)/(6.*Pi) + (Alfa*A0(MB2))/(3.*Pi) + (4*Alfa*A0(MC2))/(3.*Pi) + 
     -  (Alfa*A0(MD2))/(3.*Pi) + (Alfa*A0(ME2))/Pi + (Alfa*A0(ML2))/Pi + 
     -  (Alfa*A0(MM2))/Pi + (Alfa*A0(MS2))/(3.*Pi) + (4*Alfa*A0(MT2))/(3.*Pi) + 
     -  (4*Alfa*A0(MU2))/(3.*Pi) - (3*Alfa*A0(MW2))/(2.*Pi) + 
     -  (5*Alfa*SS*B0(SS,MW2,MW2))/(4.*Pi) - (2*Alfa*B00(SS,MB2,MB2))/(3.*Pi) - 
     -  (8*Alfa*B00(SS,MC2,MC2))/(3.*Pi) - (2*Alfa*B00(SS,MD2,MD2))/(3.*Pi) - 
     -  (2*Alfa*B00(SS,ME2,ME2))/Pi - (2*Alfa*B00(SS,ML2,ML2))/Pi - 
     -  (2*Alfa*B00(SS,MM2,MM2))/Pi - (2*Alfa*B00(SS,MS2,MS2))/(3.*Pi) - 
     -  (8*Alfa*B00(SS,MT2,MT2))/(3.*Pi) - (8*Alfa*B00(SS,MU2,MU2))/(3.*Pi) + 
     -  (3*Alfa*B00(SS,MW2,MW2))/Pi + (Alfa*SS*B1(SS,MB2,MB2))/(3.*Pi) + 
     -  (4*Alfa*SS*B1(SS,MC2,MC2))/(3.*Pi) + (Alfa*SS*B1(SS,MD2,MD2))/(3.*Pi) + 
     -  (Alfa*SS*B1(SS,ME2,ME2))/Pi + (Alfa*SS*B1(SS,ML2,ML2))/Pi + 
     -  (Alfa*SS*B1(SS,MM2,MM2))/Pi + (Alfa*SS*B1(SS,MS2,MS2))/(3.*Pi) + 
     -  (4*Alfa*SS*B1(SS,MT2,MT2))/(3.*Pi) + 
     -  (4*Alfa*SS*B1(SS,MU2,MU2))/(3.*Pi) + (Alfa*SS*B1(SS,MW2,MW2))/(2.*Pi)

      end function Agg

        double complex function AZZ(SS)
	implicit character (a-s,u-z)
	implicit double complex (t)
        double precision SS


#include "vars.h"

        AZZ= (Alfa*CW2*SS)/(6.*Pi*SW2) - (Alfa*A0(MB2))/(2.*CW2*Pi) + 
     -  (3*Alfa*A0(MB2))/(8.*CW2*Pi*SW2) + (Alfa*SW2*A0(MB2))/(3.*CW2*Pi) - 
     -  (Alfa*A0(MC2))/(CW2*Pi) + (3*Alfa*A0(MC2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SW2*A0(MC2))/(3.*CW2*Pi) - (Alfa*A0(MD2))/(2.*CW2*Pi) + 
     -  (3*Alfa*A0(MD2))/(8.*CW2*Pi*SW2) + (Alfa*SW2*A0(MD2))/(3.*CW2*Pi) - 
     -  (Alfa*A0(ME2))/(2.*CW2*Pi) + (Alfa*A0(ME2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SW2*A0(ME2))/(CW2*Pi) - (Alfa*A0(MH2))/(16.*CW2*Pi*SW2) - 
     -  (Alfa*A0(ML2))/(2.*CW2*Pi) + (Alfa*A0(ML2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SW2*A0(ML2))/(CW2*Pi) - (Alfa*A0(MM2))/(2.*CW2*Pi) + 
     -  (Alfa*A0(MM2))/(8.*CW2*Pi*SW2) + (Alfa*SW2*A0(MM2))/(CW2*Pi) - 
     -  (Alfa*A0(MS2))/(2.*CW2*Pi) + (3*Alfa*A0(MS2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SW2*A0(MS2))/(3.*CW2*Pi) - (Alfa*A0(MT2))/(CW2*Pi) + 
     -  (3*Alfa*A0(MT2))/(8.*CW2*Pi*SW2) + (4*Alfa*SW2*A0(MT2))/(3.*CW2*Pi) - 
     -  (Alfa*A0(MU2))/(CW2*Pi) + (3*Alfa*A0(MU2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SW2*A0(MU2))/(3.*CW2*Pi) + (Alfa*A0(MW2))/(4.*Pi) - 
     -  (9*Alfa*CW2*A0(MW2))/(8.*Pi*SW2) - (Alfa*SW2*A0(MW2))/(8.*CW2*Pi) - 
     -  (Alfa*A0(MZ2))/(16.*CW2*Pi*SW2) + 
     -  (3*Alfa*MB2*B0(SS,MB2,MB2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MC2*B0(SS,MC2,MC2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MD2*B0(SS,MD2,MD2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*ME2*B0(SS,ME2,ME2))/(8.*CW2*Pi*SW2) - 
     -  (Alfa*MW2*B0(SS,MH2,MZ2))/(4.*CW2**2*Pi*SW2) + 
     -  (Alfa*ML2*B0(SS,ML2,ML2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*MM2*B0(SS,MM2,MM2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MS2*B0(SS,MS2,MS2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MT2*B0(SS,MT2,MT2))/(8.*CW2*Pi*SW2) + 
     -  (3*Alfa*MU2*B0(SS,MU2,MU2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*CW2*MW2*B0(SS,MW2,MW2))/(2.*Pi*SW2) 

        AZZ = AZZ +
     -  (5*Alfa*CW2*SS*B0(SS,MW2,MW2))/(4.*Pi*SW2) - 
     -  (Alfa*MW2*SW2*B0(SS,MW2,MW2))/(2.*CW2*Pi) - 
     -  (3*Alfa*B00(SS,0D0,0D0))/(4.*CW2*Pi*SW2) +
     -  (Alfa*B00(SS,MB2,MB2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MB2,MB2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MB2,MB2))/(3.*CW2*Pi) + 
     -  (2*Alfa*B00(SS,MC2,MC2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MC2,MC2))/(4.*CW2*Pi*SW2) - 
     -  (8*Alfa*SW2*B00(SS,MC2,MC2))/(3.*CW2*Pi) + 
     -  (Alfa*B00(SS,MD2,MD2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MD2,MD2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MD2,MD2))/(3.*CW2*Pi) + 
     -  (Alfa*B00(SS,ME2,ME2))/(CW2*Pi) - 
     -  (Alfa*B00(SS,ME2,ME2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,ME2,ME2))/(CW2*Pi) + 
     -  (Alfa*B00(SS,MH2,MZ2))/(4.*CW2*Pi*SW2) + 
     -  (Alfa*B00(SS,ML2,ML2))/(CW2*Pi) - 
     -  (Alfa*B00(SS,ML2,ML2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,ML2,ML2))/(CW2*Pi) + 
     -  (Alfa*B00(SS,MM2,MM2))/(CW2*Pi) - 
     -  (Alfa*B00(SS,MM2,MM2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MM2,MM2))/(CW2*Pi) + 
     -  (Alfa*B00(SS,MS2,MS2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MS2,MS2))/(4.*CW2*Pi*SW2) - 
     -  (2*Alfa*SW2*B00(SS,MS2,MS2))/(3.*CW2*Pi) + 
     -  (2*Alfa*B00(SS,MT2,MT2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MT2,MT2))/(4.*CW2*Pi*SW2) - 
     -  (8*Alfa*SW2*B00(SS,MT2,MT2))/(3.*CW2*Pi) + 
     -  (2*Alfa*B00(SS,MU2,MU2))/(CW2*Pi) - 
     -  (3*Alfa*B00(SS,MU2,MU2))/(4.*CW2*Pi*SW2) - 
     -  (8*Alfa*SW2*B00(SS,MU2,MU2))/(3.*CW2*Pi) - 
     -  (Alfa*B00(SS,MW2,MW2))/(2.*Pi) + 
     -  (9*Alfa*CW2*B00(SS,MW2,MW2))/(4.*Pi*SW2) + 
     -  (Alfa*SW2*B00(SS,MW2,MW2))/(4.*CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,0D0,0D0))/(8.*CW2*Pi*SW2) - 
     -  (Alfa*SS*B1(SS,MB2,MB2))/(2.*CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MB2,MB2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MB2,MB2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MC2,MC2))/(CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MC2,MC2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SS*SW2*B1(SS,MC2,MC2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MD2,MD2))/(2.*CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MD2,MD2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MD2,MD2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,ME2,ME2))/(2.*CW2*Pi) + 
     -  (Alfa*SS*B1(SS,ME2,ME2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,ME2,ME2))/(CW2*Pi) - 
     -  (Alfa*SS*B1(SS,ML2,ML2))/(2.*CW2*Pi) + 
     -  (Alfa*SS*B1(SS,ML2,ML2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,ML2,ML2))/(CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MM2,MM2))/(2.*CW2*Pi) + 
     -  (Alfa*SS*B1(SS,MM2,MM2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MM2,MM2))/(CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MS2,MS2))/(2.*CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MS2,MS2))/(8.*CW2*Pi*SW2) + 
     -  (Alfa*SS*SW2*B1(SS,MS2,MS2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MT2,MT2))/(CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MT2,MT2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SS*SW2*B1(SS,MT2,MT2))/(3.*CW2*Pi) - 
     -  (Alfa*SS*B1(SS,MU2,MU2))/(CW2*Pi) + 
     -  (3*Alfa*SS*B1(SS,MU2,MU2))/(8.*CW2*Pi*SW2) + 
     -  (4*Alfa*SS*SW2*B1(SS,MU2,MU2))/(3.*CW2*Pi) + 
     -  (Alfa*CW2*SS*B1(SS,MW2,MW2))/(2.*Pi*SW2)

        end function AZZ

        double complex function AgZ(SS)
	implicit character (a-s,u-z)
	implicit double complex (t)
        double precision SS


#include "vars.h"


      AgZ= -(Alfa*CW*SS)/(6.*Pi*SW) - (Alfa*A0(MB2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SW*A0(MB2))/(6.*CW*Pi) + (Alfa*SW2*A0(MB2))/(6.*CW*Pi*SW) - 
     -  (Alfa*A0(MC2))/(2.*CW*Pi*SW) + (2*Alfa*SW*A0(MC2))/(3.*CW*Pi) + 
     -  (2*Alfa*SW2*A0(MC2))/(3.*CW*Pi*SW) - (Alfa*A0(MD2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SW*A0(MD2))/(6.*CW*Pi) + (Alfa*SW2*A0(MD2))/(6.*CW*Pi*SW) - 
     -  (Alfa*A0(ME2))/(4.*CW*Pi*SW) + (Alfa*SW*A0(ME2))/(2.*CW*Pi) + 
     -  (Alfa*SW2*A0(ME2))/(2.*CW*Pi*SW) - (Alfa*A0(ML2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SW*A0(ML2))/(2.*CW*Pi) + (Alfa*SW2*A0(ML2))/(2.*CW*Pi*SW) - 
     -  (Alfa*A0(MM2))/(4.*CW*Pi*SW) + (Alfa*SW*A0(MM2))/(2.*CW*Pi) + 
     -  (Alfa*SW2*A0(MM2))/(2.*CW*Pi*SW) - (Alfa*A0(MS2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SW*A0(MS2))/(6.*CW*Pi) + (Alfa*SW2*A0(MS2))/(6.*CW*Pi*SW) - 
     -  (Alfa*A0(MT2))/(2.*CW*Pi*SW) + (2*Alfa*SW*A0(MT2))/(3.*CW*Pi) + 
     -  (2*Alfa*SW2*A0(MT2))/(3.*CW*Pi*SW) - (Alfa*A0(MU2))/(2.*CW*Pi*SW) + 
     -  (2*Alfa*SW*A0(MU2))/(3.*CW*Pi) + (2*Alfa*SW2*A0(MU2))/(3.*CW*Pi*SW) + 
     -  (Alfa*CW*A0(MW2))/(Pi*SW) + (Alfa*CW2*A0(MW2))/(4.*CW*Pi*SW) - 
     -  (Alfa*SW2*A0(MW2))/(4.*CW*Pi*SW) - 
     -  (Alfa*CW*MW2*B0(SS,MW2,MW2))/(2.*Pi*SW) - 
     -  (5*Alfa*CW*SS*B0(SS,MW2,MW2))/(4.*Pi*SW) - 
     -  (Alfa*MW2*SW*B0(SS,MW2,MW2))/(2.*CW*Pi) + 
     -  (Alfa*B00(SS,MB2,MB2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,MB2,MB2))/(3.*CW*Pi) - 
     -  (Alfa*SW2*B00(SS,MB2,MB2))/(3.*CW*Pi*SW) + 
     -  (Alfa*B00(SS,MC2,MC2))/(CW*Pi*SW) - 
     -  (4*Alfa*SW*B00(SS,MC2,MC2))/(3.*CW*Pi) - 
     -  (4*Alfa*SW2*B00(SS,MC2,MC2))/(3.*CW*Pi*SW) + 
     -  (Alfa*B00(SS,MD2,MD2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,MD2,MD2))/(3.*CW*Pi) - 
     -  (Alfa*SW2*B00(SS,MD2,MD2))/(3.*CW*Pi*SW) + 
     -  (Alfa*B00(SS,ME2,ME2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,ME2,ME2))/(CW*Pi) - 
     -  (Alfa*SW2*B00(SS,ME2,ME2))/(CW*Pi*SW) + 
     -  (Alfa*B00(SS,ML2,ML2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,ML2,ML2))/(CW*Pi) - 
     -  (Alfa*SW2*B00(SS,ML2,ML2))/(CW*Pi*SW) + 
     -  (Alfa*B00(SS,MM2,MM2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,MM2,MM2))/(CW*Pi) - 
     -  (Alfa*SW2*B00(SS,MM2,MM2))/(CW*Pi*SW) + 
     -  (Alfa*B00(SS,MS2,MS2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SW*B00(SS,MS2,MS2))/(3.*CW*Pi) - 
     -  (Alfa*SW2*B00(SS,MS2,MS2))/(3.*CW*Pi*SW) + 
     -  (Alfa*B00(SS,MT2,MT2))/(CW*Pi*SW) - 
     -  (4*Alfa*SW*B00(SS,MT2,MT2))/(3.*CW*Pi) - 
     -  (4*Alfa*SW2*B00(SS,MT2,MT2))/(3.*CW*Pi*SW) + 
     -  (Alfa*B00(SS,MU2,MU2))/(CW*Pi*SW) - 
     -  (4*Alfa*SW*B00(SS,MU2,MU2))/(3.*CW*Pi) - 
     -  (4*Alfa*SW2*B00(SS,MU2,MU2))/(3.*CW*Pi*SW) - 
     -  (2*Alfa*CW*B00(SS,MW2,MW2))/(Pi*SW) - 
     -  (Alfa*CW2*B00(SS,MW2,MW2))/(2.*CW*Pi*SW) + 
     -  (Alfa*SW2*B00(SS,MW2,MW2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MB2,MB2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,MB2,MB2))/(6.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,MB2,MB2))/(6.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MC2,MC2))/(2.*CW*Pi*SW) + 
     -  (2*Alfa*SS*SW*B1(SS,MC2,MC2))/(3.*CW*Pi) + 
     -  (2*Alfa*SS*SW2*B1(SS,MC2,MC2))/(3.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MD2,MD2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,MD2,MD2))/(6.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,MD2,MD2))/(6.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,ME2,ME2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,ME2,ME2))/(2.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,ME2,ME2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,ML2,ML2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,ML2,ML2))/(2.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,ML2,ML2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MM2,MM2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,MM2,MM2))/(2.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,MM2,MM2))/(2.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MS2,MS2))/(4.*CW*Pi*SW) + 
     -  (Alfa*SS*SW*B1(SS,MS2,MS2))/(6.*CW*Pi) + 
     -  (Alfa*SS*SW2*B1(SS,MS2,MS2))/(6.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MT2,MT2))/(2.*CW*Pi*SW) + 
     -  (2*Alfa*SS*SW*B1(SS,MT2,MT2))/(3.*CW*Pi) + 
     -  (2*Alfa*SS*SW2*B1(SS,MT2,MT2))/(3.*CW*Pi*SW) - 
     -  (Alfa*SS*B1(SS,MU2,MU2))/(2.*CW*Pi*SW) + 
     -  (2*Alfa*SS*SW*B1(SS,MU2,MU2))/(3.*CW*Pi) + 
     -  (2*Alfa*SS*SW2*B1(SS,MU2,MU2))/(3.*CW*Pi*SW) - 
     -  (Alfa*CW*SS*B1(SS,MW2,MW2))/(2.*Pi*SW)

      end function AgZ
