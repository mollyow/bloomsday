# Asymptotic Statistics Quarto Book

This directory contains a Quarto book of study notes for selected chapters of A. W. van der Vaart's *Asymptotic Statistics*.

## Requirements

- [Quarto](https://quarto.org/docs/get-started/)

Check the installation with:

```sh
quarto --version
```

## Preview locally

From this directory, run:

```sh
quarto preview
```

Quarto will start a local server and rebuild changed pages automatically.

## Render the book

```sh
quarto render
```

The rendered site is written to `_book/`, beginning at `_book/index.html`. The output directory and Quarto's local cache are excluded from Git.

## Project organization

- `_quarto.yml` defines the book and chapter order.
- `index.qmd` contains the reading plan.
- `chapter-*.qmd` contains the chapter notes.
- `glossary.qmd` collects shared notation.
- `problems.qmd` contains worked problems.
- `images/` contains book assets referenced by the notes.

The source book's chapter and result numbers are written explicitly because this project includes a nonconsecutive selection of chapters. Major cross-references use stable anchors so that links survive nearby heading edits.

## Local source material

The original Notion export and reference PDFs are retained locally as source material but are intentionally excluded from Git. The Quarto files in this directory are the maintained version of the notes.
