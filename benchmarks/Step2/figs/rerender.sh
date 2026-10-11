#!/bin/bash
# gnuplot -> fig -> pdf/eps in figs/ (as make_examples.render): rerender.sh stem
set -e
cd "$(dirname "$0")"
stem=$1
gnuplot $stem.gp
figtex2pdf.sh $stem.fig > /dev/null
mv $stem.pdf ../$stem.pdf
pdftops -eps ../$stem.pdf ../$stem.eps
echo "figs/$stem.pdf, .eps"
