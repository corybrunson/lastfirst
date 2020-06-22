# Fixed and adaptive landmark sets for finite metric spaces

This folder contains the base manuscript, which consists in the following files:

* `fixed-adaptive.md`: The Markdown file to edit directly.
* `fixed-adaptive.bib`: The Mendeley-generated bibliography; this file should not be edited directly.
* `format.sty`: A file containing LaTeX formatting used to render `fixed-adaptive.md`

As commented at the end of `fixed-adaptive.md`, the LaTeX file and PDF can be generated using the following Pandoc commands:

```sh
# To generate the LaTeX file, execute the following:
pandoc fixed-adaptive.md \
-o fixed-adaptive.tex \
-s \
--number-sections \
--bibliography=fixed-adaptive.bib
# To generate the PDF directly, execute the following:
pandoc fixed-adaptive.md \
-o fixed-adaptive.pdf \
-s \
--number-sections \
--bibliography=fixed-adaptive.bib
```

The Pandoc-generated files `fixed-adaptive.tex` and `fixed-adaptive.pdf` are ignored in order to avoid unnecessary conflicts.
