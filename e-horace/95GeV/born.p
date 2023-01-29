set title "Z boson invariant mass distribution"
set xlabel "Invariant mass [Gev]"
set ylabel "Counts"

plot "alpha/tev_z_mu_minv_oal_norec_200.dat" title "ew=0 alpha" with errorbars, \
     "testalpha/tev_z_mu_minv_oal_norec_200.dat" title "ew=1 alpha" with errorbars, \
     "born/tev_z_mu_minv_born_norec_200.dat" title "ew=0 born" with errorbars, \
     "testborn/tev_z_mu_minv_born_norec_200.dat" title "ew=1 born" with errorbars, \
     "testoalpha/1/tev_z_mu_minv_oal_norec_200.dat" title "ew=1 old" with errorbars, \
     "testoalpha/0/tev_z_mu_minv_oal_norec_200.dat" title "ew=0 old" with errorbars


