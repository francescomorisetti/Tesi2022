
set title "Energia nel centro di massa=95GeV"
set xlabel "Massa invariante [Gev]"

set key top center
set xrange[84:95]

plot "born/tev_z_mu_minv_born_norec_200.dat" title "born" with errorbars lw 0.5 lt rgb 'black', \
     "alpha/tev_z_mu_minv_oal_norec_200.dat" title "alpha" with errorbars lt rgb 'red', \
     "exp/tev_z_mu_minv_best_norec_200.dat" title "exp" with errorbars lt rgb 'dark-green', \
     "oldexp/tev_z_mu_minv_best_norec_200.dat" title "old-exp" with errorbars lt rgb 'blue', \
     "oldalpha/tev_z_mu_minv_oal_norec_200.dat" title "old-alpha" with errorbars lw 0.2 lt rgb 'gold'
     


