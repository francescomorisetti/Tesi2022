
set title "Z boson invariant mass distribution"
set xlabel "Invariant mass [Gev]"
set ylabel "Counts"

plot "100GeV/tev_z_mu_minv_born_norec_200.dat" title "Q=100GeV" lc rgb 'red' w errorbars, \
     "50GeV/tev_z_mu_minv_born_norec_200.dat" title "Q=50GeV" lc rgb 'blue' w errorbars, \
     "25GeV/tev_z_mu_minv_born_norec_200.dat" title "Q=25GeV" lc rgb 'green' w errorbars   
      

