set output "surf_ac.fig"
set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set linetype 1 lw 2 pt 7 ps 1.1 lc rgb "#000000"
set linetype 2 lw 2 pt 5 ps 1.1 lc rgb "#E69F00"
set linetype 3 lw 2 pt 9 ps 1.1 lc rgb "#56B4E9"
set linetype 4 lw 2 pt 13 ps 1.1 lc rgb "#009E73"
set linetype 5 lw 2 pt 15 ps 1.1 lc rgb "#D55E00"
set linetype 6 lw 2 pt 11 ps 1.1 lc rgb "#0072B2"
set linetype 7 lw 2 pt 3 ps 1.1 lc rgb "#CC79A7"
set xlabel "$a/t$"
set ylabel "$a/c$"
set key top left Left reverse nobox spacing 1.4
plot "surf_ac.txt" u 5:7 w l lt 5 t "Newman--Raju", "surf_ac.txt" u 5:4 w p lt 1 t "computed"
unset output
