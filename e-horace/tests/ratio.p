set title "EW scheme comparison"
set xlabel "Invariant mass [GeV]"
set ylabel "Ratio"

plot '<paste 1/tev_z_mu_minv_oal_norec_200.dat 0/tev_z_mu_minv_oal_norec_200.dat' using 1:($5/$2) w l
