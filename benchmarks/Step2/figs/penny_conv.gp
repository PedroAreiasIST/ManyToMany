set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set output "penny_conv.fig"
set linetype 1 lw 2 pt 7 ps 1.1 lc rgb "#000000"
set linetype 2 lw 2 pt 5 ps 1.1 lc rgb "#E69F00"
set linetype 3 lw 2 pt 9 ps 1.1 lc rgb "#56B4E9"
set linetype 4 lw 2 pt 13 ps 1.1 lc rgb "#009E73"
set linetype 5 lw 2 pt 15 ps 1.1 lc rgb "#D55E00"
set linetype 6 lw 2 pt 11 ps 1.1 lc rgb "#0072B2"
set linetype 7 lw 2 pt 3 ps 1.1 lc rgb "#CC79A7"
set logscale xy
set xlabel "$h_{\\min}/a$"
set ylabel "$\\max_\\theta |K-K_{\\mathrm{ref}}|/K_0$ [\\%]"
set xrange [0.025:0.15]
set yrange [0.1:5]
set key outside right top Left reverse nobox spacing 1.4
plot "penny_conv.txt" u 1:2 w lp lt 1 t "$\\beta=0$: $K_{\\mathrm{I}}$", \
     "" u 1:3 w lp lt 2 t "$\\beta=30^\\circ$: $K_{\\mathrm{I}}$", \
     "" u 1:4 w lp lt 5 t "$\\beta=30^\\circ$: $K_{\\mathrm{II}}$", \
     "" u 1:5 w lp lt 3 t "$\\beta=30^\\circ$: $K_{\\mathrm{III}}$", \
     12*x w l lt 7 dt 2 t "slope 1", 80*x**2 w l lt 7 dt 3 t "slope 2"
unset output
