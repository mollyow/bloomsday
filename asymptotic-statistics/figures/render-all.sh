#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

mkdir -p images/generated figures/build
Rscript figures/R/render-figures.R

for source in figures/tikz/*.tex; do
  name="$(basename "$source" .tex)"
  pdflatex -interaction=batchmode -halt-on-error \
    -output-directory figures/build "$source" >/dev/null
  pdftocairo -svg "figures/build/${name}.pdf" \
    "images/generated/${name}.svg"
done

rm -f figures/build/*.aux figures/build/*.dvi figures/build/*.log \
  figures/build/*.pdf figures/build/*.xdv

echo "Rendered figures to images/generated"
