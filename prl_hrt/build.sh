#!/bin/bash
# Regenerate the figures and build letter.pdf and supplement.pdf.
set -e
cd "$(dirname "$0")"
( cd fig
  python3 fig1_wormhole.py fig1_wormhole.pdf
  python3 fig2_mechanism.py fig2_mechanism.pdf
  python3 fig3_numerics.py ../../fullcone fig3_numerics.pdf )
for f in letter supplement; do
  pdflatex -interaction=nonstopmode $f.tex >/dev/null
  bibtex $f >/dev/null || true
  pdflatex -interaction=nonstopmode $f.tex >/dev/null
  bibtex $f >/dev/null || true   # second pass: citations inside the SM note of the bibliography
  pdflatex -interaction=nonstopmode $f.tex >/dev/null
  pdflatex -interaction=nonstopmode $f.tex >/dev/null
  grep -q 'undefined' $f.log && { echo "undefined references in $f"; exit 1; } || true
  grep -q 'Overfull' $f.log && echo "warning: overfull boxes in $f" || true
done
rm -f *.aux *.log *.blg *.out *Notes.bib fig/__pycache__ -r
echo "built letter.pdf ($(pdfinfo letter.pdf | awk '/Pages/{print $2}') pp) and supplement.pdf ($(pdfinfo supplement.pdf | awk '/Pages/{print $2}') pp)"
