
set title "Z boson invariant mass distribution"
set xlabel "Invariant mass [Gev]"
set ylabel "Counts"

plot "born/tev_z_mu_minv_born_norec_200.dat" title "born" with errorbars, \
     "alpha/tev_z_mu_minv_oal_norec_200.dat" title "alpha" with errorbars, \
     "exp/tev_z_mu_minv_best_norec_200.dat" title "expo" with errorbars, \
     "oldexp/tev_z_mu_minv_best_norec_200.dat" title "old-expo" with errorbars, \
     "oldalpha/tev_z_mu_minv_oal_norec_200.dat" title "old-alpha" with errorbars
     


