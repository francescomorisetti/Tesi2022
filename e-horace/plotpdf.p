set title "NLL DELTA MSBAR, Q=100GeV, {/Symbol g}=1"
set xlabel "x"
set ylabel "f(x)"

set log y

set key top center

plot "elpdf.dat" w l title "{e^-}" lt rgb 'red' lw 1.5, \
     "pospdf.dat" w l title "{e^+}" lt rgb 'forest-green' lw 1.5, \
     "phpdf.dat" w l title "{/Symbol g}" lt rgb 'blue' lw 1.5
