set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set output "penny_K.fig"
set linetype 1 lw 2 pt 7 ps 1.1 lc rgb "#000000"
set linetype 2 lw 2 pt 5 ps 1.1 lc rgb "#E69F00"
set linetype 3 lw 2 pt 9 ps 1.1 lc rgb "#56B4E9"
set linetype 4 lw 2 pt 13 ps 1.1 lc rgb "#009E73"
set linetype 5 lw 2 pt 15 ps 1.1 lc rgb "#D55E00"
set linetype 6 lw 2 pt 11 ps 1.1 lc rgb "#0072B2"
set linetype 7 lw 2 pt 3 ps 1.1 lc rgb "#CC79A7"
set xlabel "Angle from the shear direction $\\theta$ [deg]"
set ylabel "$K/K_0$"
set xrange [0:360]
set xtics 0,90,360
set yrange [-0.75:1.0]
set key outside right top Left reverse nobox spacing 1.4
plot "penny_Kref.txt" u 1:2 w l lt 1 dt 1 t "$K_{\\mathrm{I}}$, Kassir--Sih", \
     "penny_K.txt" u 1:2 w p lt 1 t "$K_{\\mathrm{I}}$, computed", \
     "penny_Kref.txt" u 1:3 w l lt 5 t "$K_{\\mathrm{II}}$, Kassir--Sih", \
     "penny_K.txt" u 1:3 w p lt 5 t "$K_{\\mathrm{II}}$, computed", \
     "penny_Kref.txt" u 1:4 w l lt 3 t "$K_{\\mathrm{III}}$, Kassir--Sih", \
     "penny_K.txt" u 1:4 w p lt 3 t "$K_{\\mathrm{III}}$, computed"
unset output
