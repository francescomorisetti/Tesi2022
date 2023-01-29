
set title "Z boson invariant mass distribution"
set xlabel "Invariant mass [Gev]"
set ylabel "Counts"

plot "alpha/11/tev_z_mu_minv_oal_norec_200.dat" title "alpha ew=0" with errorbars, \
     "test/alpha/tev_z_mu_minv_oal_norec_200.dat" title "alpha ew=1" with errorbars, \
     "born/11/tev_z_mu_minv_born_norec_200.dat" title "born ew=0" with errorbars, \
     "test/born/tev_z_mu_minv_born_norec_200.dat" title "born ew=1" with errorbars
     
