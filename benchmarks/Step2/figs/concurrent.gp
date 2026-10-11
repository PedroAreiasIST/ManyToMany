# Readers during a crack sequence on 10^6 elements.
set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set output "concurrent.fig"
set linetype 1 lw 2 pt 7  ps 1.2 lc rgb "#000000"
set linetype 6 lw 2 pt 15 ps 1.2 lc rgb "#D55E00"
set logscale y
set format y "$10^{%T}$"
set xlabel "Reader threads"
set ylabel "Reads per second (all readers)"
set xrange [0.5:16]
set key inside bottom right Right nobox spacing 1.4
plot "concurrent_step2.txt" u 1:2 w lp lt 1 t "snapshots (lock-free)", \
     "concurrent_mm3.txt" u 1:2 w lp lt 6 t "MM3 reader--writer lock"
unset output
