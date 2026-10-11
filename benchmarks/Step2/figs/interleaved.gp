# Recursive bisection, queries interleaved with edits: cost per bisection.
set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set output "interleaved.fig"
set linetype 1 lw 2 pt 7  ps 1.2 lc rgb "#000000"
set linetype 6 lw 2 pt 15 ps 1.2 lc rgb "#D55E00"
set logscale xy
set format x "$10^{%T}$"
set format y "$10^{%T}$"
set xlabel "Elements $n_E$ (initial mesh)"
set ylabel "Time per bisection [$\\mu$s]"
set key inside top left Left reverse nobox spacing 1.4
plot "interleaved.txt" u 1:4 w lp lt 6 t "MM3 (recompute on query)", \
     "" u 1:3 w lp lt 1 t "incremental"
unset output
