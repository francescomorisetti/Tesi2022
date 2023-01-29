set title "Invariant mass ratio"
set xlabel "Invariant mass [GeV]"
set ylabel "Ratio"

set xrange [88:94]


plot '<paste oldalpha/tev_z_mu_minv_oal_norec_200.dat born/tev_z_mu_minv_born_norec_200.dat' using 1:($2/$5) w l title "Oldalpha/Born" lt rgb 'red' lw 1.5, \
     '<paste oldexp/tev_z_mu_minv_best_norec_200.dat born/tev_z_mu_minv_born_norec_200.dat' using 1:($2/$5) w l title "Oldexp/Born" lt rgb 'green'lw 1.5, \
     '<paste alpha/tev_z_mu_minv_oal_norec_200.dat born/tev_z_mu_minv_born_norec_200.dat' using 1:($2/$5) w l title "Alpha/Born" lt rgb 'orange' lw 1.5, \
     '<paste exp/tev_z_mu_minv_best_norec_200.dat born/tev_z_mu_minv_born_norec_200.dat' using 1:($2/$5) w l title "Expo/Born" lt rgb 'blue' lw 1.5, \
     '<paste born/tev_z_mu_minv_born_norec_200.dat born/tev_z_mu_minv_born_norec_200.dat' using 1:($2/$5) w l title "Born/Born" lt rgb 'black' lw 1.5
