
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flod

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/flod)](https://CRAN.R-project.org/package=flod)
[![Codecov test
coverage](https://codecov.io/gh/hugoakerstrand/flod/graph/badge.svg)](https://app.codecov.io/gh/hugoakerstrand/flod)
[![R-CMD-check](https://github.com/hugoakerstrand/flod/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hugoakerstrand/flod/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Introduction

flod is an R package for statistical learning on flow cytometry data in
R. Its primary concern is getting your flow cytometry data accessible to
R’s rich ecosystem for modeling, for example
[tidymodels](https://www.tidymodels.org/).

Additionally, flod also provides a few convenience functions for feature
engineering, visualizing, and reporting of flow cytometry data.

## Background

Flow cytometry quantifies light scattering properties of a cell to
measure its size, granularity, and any biomolecules that have been
pre-labelled by the experimentalist.

These values are stored in the .fcs file format that contains raw
expression values of single events, sample meta data, and run
information.

The data is then analyzed by sequentially filtering down expression and
density of events in one or two dimensions at a time, typically using
proprietary point-and-click software.

The power of flow cytometry is its reasonable performance metrics across
the board: the operator can analyze thousands of cells per second for
several markers of various nature. Combined with a relatively low
cost-per-run and ease of operation, flow cytometry has become a pillar
of modern cell based research.

## Package motivation

flod makes flow cytometry data accessible to statistical analysis and
modelling, by functions to read in .fcs files and the feature
engineering expected by flow cytometry experts.

More specifically, it was developed to adress the following challenges
of standard flow cytometry data analysis: - **low reproducibility**:
while ease of use software has made flow cytometry data analysis
accessible, it has a lot of arbirtary choices by the operator that makes
reproducibility low. - **time consuming**: data analysis is time
consuming and routine analysis becomes a big time sink. - **unsuited for
production environment**: analysis is in proprietary software is not
approved by regulatory agencies such as FDA, as of 2026-04-30. -
**inaccessible to more advanced statistical tools**: lifting data
analysis into R puts it in lockstep with current state-of-art tools for
statistical modelling.
