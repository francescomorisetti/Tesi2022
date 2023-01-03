
set title "Z boson invariant mass distribution"
set xlabel "Invariant mass [Gev]"
set ylabel "Counts"

plot "0/tev_z_mu_minv_oal_norec_200.dat" title "Pure alpha" with errorbars, \
     "1/tev_z_mu_minv_oal_norec_200.dat" title "G_\mu" with errorbars


