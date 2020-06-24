---
title: "Fixed and adaptive landmark sets for finite metric spaces"
author:
  - Jason Cory Brunson
  - Yara Skaf
output:
  pdf_document:
    fig_caption: yes
    keep_tex: yes
    number_sections: no
header-includes:
  - \usepackage{natbib}
  - \usepackage{format}
  - \setcitestyle{numbers}
---

<!--
header-includes:
  - \usepackage{lineno}
  - \linenumbers
  - \doublespacing
-->

\pagebreak

# Introduction

## Notation and Conventions

$(X, d_X)$ will refer to a finite pseudometric space with point set $X$ and pseudometric $d_X:X\times X\to\R_{\geq 0}$; $(X,d_X)$ may be shortened to $X$, and $d_X$ to $d$, when clear from context. The cardinality of $X$ is denoted $\abs{X}$.
$f:X \to Y$ will denote a morphism of pseudometric spaces, i.e. a morphism of sets satisfying $d_X(x,y)\geq d_Y(f(x),f(y))$.
If $x\neq y$ but $d(x,y)=0$ then $x$ and $y$ are said to be co-located.
If $d(x,y)=d(x,z)\implies y=z$, then $X$ is considered to be in general position, even if $d(x,y)=d(z,w)\nimplies \{x,y\}=\{z,w\}$. Note that this is condition implies that $X$ is Hausdorff ($d(x,y)=0\implies x=y$).

We use the open ball notation $B_{\eps}(x)$ for the set of points less than distance $\eps$ from a central point $x$; that is, $B_{\eps}(x) = \{ y \mid d(x,y) < \eps \}$.
We use an overline to include points exactly distance $\eps$ from $x$: $\cl{B_{\eps}}(x) = \{ y \mid d(x,y) \leq \eps \}$.
When $\abs{B_\eps(x)}\geq k$ and $\eps'<\eps\implies\abs{B_\eps(x)}<k$, then we say that $N^+_k(x)=B_\eps(x)$ is the $k$-nearest neighborhood of $x$. When $X$ is in general position, $\abs{N_k(x)}=k$.

For convenience, we take the natural numbers $\N$ to include $0$.

## Background

Motivations for the maxmin procedure

* vertex set for a witness complex construction [@DeSilva2004]
    * simplices defined inductively from nearest neighborhoods of points
    * heuristically optimal in-neighborhood cover, locally well-separated
    * avoid density bias that might accentuate accidental features
* cardinality reduction of a point cloud in preparation for mapper [@Singh2007]
* optimal fixed-radius ball cover for a point cloud [@Dlotko2019]
    * nerve that approximates mapper
    * filtration that approximate persistent homology
* less geometric, more topological cover for mapper constructions
    * relies on meaningful distance metric

## Motivation

Motivations for the lastfirst procedure

* optimal fixed-cardinality out-neighborhood cover for a point cloud
    * does not subsume nearby points with high multiplicity
    * accommodates relative distance or similarity measures
    * easily coupled with nearest neighbor–based predictive models
* basis for 2-parameter counterpart (out- versus in-) to witness complexes?

# Procedures

## Maxmin procedure

## Lastfirst procedure

### Rank distances

### Tie handling

# Evaluations

## Data sets

### Simulated data

### Empirical data

## Benchmark tests

## Robustness analysis

## Stability analysis

# Applications

## Data visualization

### CCU data

### CTD data

## Interpolative prediction

### CCU data

### CTD data

# References

<!--
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
-->
