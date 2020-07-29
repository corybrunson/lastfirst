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

## Background

## Motivation

# Procedures

## Maxmin procedure

## Lastfirst procedure

### Rank distances

### Tie handling

# Evaluations

## Data sets

### Simulated data

_Choose one (at least moderately) large data set for each combination of the following properties: high- versus low-dimensional, with and without high multiplicities. Have the simulated data sets exceed/bookend the empirical data sets on both sides of each property range._

### Empirical data

* [MIMIC-III Cardiac Care Unit](https://mimic.physionet.org/mimictables/transfers/) under one or both similarity measures
* [México COVID-19 Comunicado Técnico Diario](https://www.gob.mx/salud/documentos/coronavirus-covid-19-comunicado-tecnico-diario-238449) under one similarity measure

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
--bibliography=../lastfirst.bib
# To generate the PDF directly, execute the following:
pandoc fixed-adaptive.md \
-o fixed-adaptive.pdf \
-s \
--number-sections \
--bibliography=../lastfirst.bib
-->
