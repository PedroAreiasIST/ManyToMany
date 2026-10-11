#!/bin/bash
# Lightweight corner-crack figure: each panel (seam<N>.pdf from strip.py + seam.py + epstopdf) is rendered at 600 dpi,
# cropped (trim 28 23 68 24 pt) and stored as JPEG (quality 93); cs_jpg.tex composes the 2x2 PDF (about 0.4 MB instead of 11 MB).
set -e
for s in 5 8 10 12; do
  pdftoppm -r 304 -x $((28*304/72)) -y $((24*304/72)) -W $((444*304/72)) -H $((488*304/72)) -png seam$s.pdf p$s
  python3 -c "from PIL import Image,glob; import glob as g; Image.open(g.glob('p$s-*.png')[0]).convert('RGB').save('p$s.jpg',quality=93,subsampling=0,optimize=True)" 2>/dev/null || \
  python3 -c "from PIL import Image; import glob; Image.open(glob.glob('p$s-*.png')[0]).convert('RGB').save('p$s.jpg',quality=93,subsampling=0,optimize=True)"
done
pdflatex -interaction=nonstopmode cs_jpg.tex
