set title "Z boson rapidity distribution"
set xlabel "Rapidity [Gev]"
set ylabel "Counts"


plot "100GeV/tev_z_mu_yv_born_norec_200.dat" title "Q=100GeV" lc rgb 'red' w errorbars, \
     "50GeV/tev_z_mu_yv_born_norec_200.dat" title "Q=50GeV" lc rgb 'blue' w errorbars, \
     "25GeV/tev_z_mu_yv_born_norec_200.dat" title "Q=25GeV" lc rgb 'green' w errorbars
