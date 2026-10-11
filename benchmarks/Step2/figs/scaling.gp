# Per-step cost of the crack advance against mesh size: incremental maintenance vs recomputation.
set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set output "scaling.fig"
set linetype 1 lw 2 pt 7  ps 1.2 lc rgb "#000000"
set linetype 2 lw 2 pt 9  ps 1.2 lc rgb "#E69F00"
set linetype 3 lw 2 pt 5  ps 1.2 lc rgb "#56B4E9"
set linetype 4 lw 2 pt 13 ps 1.2 lc rgb "#009E73"
set linetype 5 lw 2 pt 11 ps 1.2 lc rgb "#0072B2"
set linetype 6 lw 2 pt 15 ps 1.2 lc rgb "#D55E00"
set logscale xy
set format x "$10^{%T}$"
set format y "$10^{%T}$"
set xlabel "Elements $n_E$"
set ylabel "Time per crack step [ms]"
set key inside left top Left reverse nobox spacing 1.3
set xrange [6e3:1.6e7]
plot "scaling_crack.txt" u 1:6 w lp lt 6 t "MM3 recomputation", \
     "" u 1:8 w lp lt 3 t "tuned recomputation, all (16 thr.)", \
     "" u 1:7 w lp lt 4 t "tuned recomputation, MM3 set (16 thr.)", \
     "" u 1:4 w lp lt 1 t "incremental, persistent", \
     "" u 1:5 w lp lt 2 t "incremental, transient"
unset output
