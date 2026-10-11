set terminal fig color big textspecial fontsize 14 linewidth 2.0 depth 10
set output "penny_cod.fig"
set linetype 1 lw 2 pt 7 ps 1.1 lc rgb "#000000"
set linetype 2 lw 2 pt 5 ps 1.1 lc rgb "#E69F00"
set linetype 3 lw 2 pt 9 ps 1.1 lc rgb "#56B4E9"
set linetype 4 lw 2 pt 13 ps 1.1 lc rgb "#009E73"
set linetype 5 lw 2 pt 15 ps 1.1 lc rgb "#D55E00"
set linetype 6 lw 2 pt 11 ps 1.1 lc rgb "#0072B2"
set linetype 7 lw 2 pt 3 ps 1.1 lc rgb "#CC79A7"
set xlabel "$r/a$"
set ylabel "Displacement jump$/(8(1-\\nu^2)\\sigma a/(\\pi E))$"
set xrange [0:1]
set yrange [0:0.8]
set key outside right top Left reverse nobox spacing 1.4
sn = 0.7500000000
st = 0.5094267081
plot "penny_cod.txt" u 1:2 w p lt 1 pt 7 ps 0.4 t "opening, computed", \
     sn*sqrt(1-x**2) w l lt 1 t "opening, Kassir--Sih", \
     "penny_cod.txt" u 1:3 w p lt 5 pt 7 ps 0.4 t "sliding, computed", \
     st*sqrt(1-x**2) w l lt 5 t "sliding, Kassir--Sih"
unset output
