set output "surf_fronts.fig"
set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set linetype 1 lw 2 pt 7 ps 1.1 lc rgb "#000000"
set linetype 2 lw 2 pt 5 ps 1.1 lc rgb "#E69F00"
set linetype 3 lw 2 pt 9 ps 1.1 lc rgb "#56B4E9"
set linetype 4 lw 2 pt 13 ps 1.1 lc rgb "#009E73"
set linetype 5 lw 2 pt 15 ps 1.1 lc rgb "#D55E00"
set linetype 6 lw 2 pt 11 ps 1.1 lc rgb "#0072B2"
set linetype 7 lw 2 pt 3 ps 1.1 lc rgb "#CC79A7"
set size ratio -1
set xlabel "$x$"
set ylabel "depth $y+t$"
set xtics 0.1
set ytics 0.1
set yrange [0:*] reverse
set key off
plot for [i=1:12] "surf_fronts.txt" index i u 1:2 w l lt 3 lw 1.5, "surf_fronts.txt" index 0 u 1:2 w l lt 1 lw 2.5
unset output
