set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set output "growth_K.fig"
set linetype 1 lw 2 pt 7 ps 1.1 lc rgb "#000000"
set linetype 2 lw 2 pt 5 ps 1.1 lc rgb "#E69F00"
set linetype 3 lw 2 pt 9 ps 1.1 lc rgb "#56B4E9"
set linetype 4 lw 2 pt 13 ps 1.1 lc rgb "#009E73"
set linetype 5 lw 2 pt 15 ps 1.1 lc rgb "#D55E00"
set linetype 6 lw 2 pt 11 ps 1.1 lc rgb "#0072B2"
set linetype 7 lw 2 pt 3 ps 1.1 lc rgb "#CC79A7"
set xlabel "Angle along the front $\\psi$ [deg]"
set ylabel "$K_{\\mathrm{I}}/(\\sigma\\sqrt{\\pi b_0})$"
set xrange [0:360]
set xtics 0,90,360
set yrange [0.6:1.2]
set key inside top center horizontal maxcols 2 Left reverse nobox spacing 1.2 width -1
plot "growth_K_0.txt" u 1:2 w lp lt 2 t "step 0", \
     "growth_K_4.txt" u 1:2 w lp lt 3 t "step 4", \
     "growth_K_8.txt" u 1:2 w lp lt 4 t "step 8", \
     "growth_K_0.txt" u 1:3 w l lt 1 t "Irwin, step 0"
unset output
