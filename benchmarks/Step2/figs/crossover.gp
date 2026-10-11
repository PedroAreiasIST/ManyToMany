# Batch of random node duplications on 10^6 elements: incremental vs recomputation.
set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set output "crossover.fig"
set linetype 1 lw 2 pt 7  ps 1.2 lc rgb "#000000"
set linetype 3 lw 2 pt 5  ps 1.2 lc rgb "#56B4E9"
set linetype 4 lw 2 pt 13 ps 1.2 lc rgb "#009E73"
set linetype 6 lw 2 pt 15 ps 1.2 lc rgb "#D55E00"
set logscale xy
set format x "$10^{%T}$"
set format y "$10^{%T}$"
set xrange [0.7:1.5e6]
set xlabel "Edits in the batch"
set ylabel "Time [ms]"
set key inside left top Left reverse nobox spacing 1.3
plot "crossover.txt" u 1:3 w lp lt 1 t "incremental, persistent", \
     "" u 1:4 w l lt 6 t "MM3 recomputation", \
     "" u 1:7 w l lt 3 t "tuned recomputation, all (16 thr.)", \
     "" u 1:6 w l lt 3 dt 2 t "tuned recomputation, all (1 thr.)", \
     "" u 1:5 w l lt 4 t "tuned recomputation, MM3 set (16 thr.)"
unset output
