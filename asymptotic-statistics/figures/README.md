# Reproducible figures

Run from the book directory:

```sh
bash figures/render-all.sh
```

The R sources create statistical plots with fixed seeds through `svglite`.
Standalone TikZ sources are compiled with `pdflatex` and converted with
Poppler's `pdftocairo`. The committed SVG files in `images/generated/` let `quarto render`
remain independent of R and LaTeX execution.

The shared visual language is blue for target or parametric objects, orange for
nuisance or alternative objects, green for efficient residuals, and gray for
approximations or supporting structure.

For full-resolution PNG previews during visual QA, run the R source with an
explicit temporary output directory, for example:

```sh
FIGURE_QA_DIR=/tmp/bloomsday-figure-qa Rscript figures/R/render-figures.R
```

This optional preview pass does not add PNGs to the book.
