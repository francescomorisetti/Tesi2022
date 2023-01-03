set title "Invariant mass ratio"
set xlabel "Invariant mass [GeV]"
set ylabel "Ratio"

plot '<paste oldalpha/tev_z_mu_minv_oal_norec_200.dat born/tev_z_mu_minv_born_norec_200.dat' using 1:($5/$2) w l title "Oldalpha/Born", \
     '<paste oldexp/tev_z_mu_minv_best_norec_200.dat born/tev_z_mu_minv_born_norec_200.dat' using 1:($5/$2) w l title "Oldexp/Born", \
     '<paste alpha/tev_z_mu_minv_oal_norec_200.dat born/tev_z_mu_minv_born_norec_200.dat' using 1:($5/$2) w l title "Alpha/Born", \
     '<paste expo/tev_z_mu_minv_best_norec_200.dat born/tev_z_mu_minv_born_norec_200.dat' using 1:($5/$2) w l title "Expo/Born", \
     '<paste born/tev_z_mu_minv_born_norec_200.dat born/tev_z_mu_minv_born_norec_200.dat' using 1:($5/$2) w l title "Born/Born"
