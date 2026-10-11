set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set output "cutpar.fig"
set linetype 1 lw 2 pt 7 ps 1.1 lc rgb "#000000"
set linetype 2 lw 2 pt 5 ps 1.1 lc rgb "#E69F00"
set linetype 3 lw 2 pt 9 ps 1.1 lc rgb "#56B4E9"
set linetype 4 lw 2 pt 13 ps 1.1 lc rgb "#009E73"
set linetype 5 lw 2 pt 15 ps 1.1 lc rgb "#D55E00"
set linetype 6 lw 2 pt 11 ps 1.1 lc rgb "#0072B2"
set linetype 7 lw 2 pt 3 ps 1.1 lc rgb "#CC79A7"
set logscale x 2
set xlabel "Threads"
set ylabel "Speed-up over the sequential cut"
set xrange [0.8:160]
set yrange [0:14]
set key inside left top Left reverse nobox spacing 1.3
plot "cutpar_0.txt" u 1:2 w lp lt 1 t "9950X, $10^7$ el., 64 groups", \
     "cutpar_1.txt" u 1:2 w lp lt 2 t "9950X, $10^7$ el., 512 groups", \
     "cutpar_2.txt" u 1:2 w lp lt 3 t "EPYC 7742, $10^7$ el., 64 groups", \
     "cutpar_3.txt" u 1:2 w lp lt 4 t "EPYC 7742, $10^7$ el., 512 groups", \
     "cutpar_4.txt" u 1:2 w lp lt 5 t "EPYC 7742, $2.5\\times10^7$ el., 512 groups"
unset output
