# CloudClean.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/andrew-saydjari/CloudClean.jl)
[![Build Status](https://github.com/andrew-saydjari/CloudClean.jl/workflows/Unit%20test/badge.svg)](https://github.com/andrew-saydjari/CloudClean.jl/actions)
[![Coverage Status](https://codecov.io/github/andrew-saydjari/CloudClean.jl/coverage.svg?branch=main)](https://codecov.io/github/andrew-saydjari/CloudClean.jl?branch=main)

Image infilling algorithm with focus on statistically accuracy

## Installation

**CloudClean** is a registered package so a stable version can be installed using `Pkg.add`.

```julia
import Pkg
Pkg.add("CloudClean")
```

For the most recent development version, install directly from the GitHub

```julia
import Pkg
Pkg.add(url="https://github.com/andrew-saydjari/CloudClean.jl")
```

## Description

By leveraging the local pixel-pixel covarariance structure in an image, CloudClean attempts to predict the values of missing pixels based on near-by unmasked pixels. The user needs only to provide an image and a mask of "bad" pixels and the choice of a few hyper parameters. CloudClean has two main operating modes:

1. proc_discrete, which attempts to fill in masked data in a subregion centered on a discrete list of input points
2. proc_continuous, which infills arbitrarily shaped and centered masks. The latter is like "Photoshop" for images with correlated structure.

This code is based heavily upon (and is in some sense a simplification of) [CloudCovErr.jl](https://github.com/andrew-saydjari/CloudCovErr.jl).

## Example

An example of the quality of the prediction for missing pixels is demonstrated on this image from the WISE 12 μm dust map. More examples, notebooks, and documentation are in process.

[!["WISE infill example"][infill-img]][infill-url]

## Contributing and Questions

This is a new piece of software. [Filing an
issue](https://github.com/andrew-saydjari/CloudClean.jl/issues/new) to report a
bug or request a feature is extremely valuable in helping us prioritize what to work on, so don't hesitate.

## Table of Contents

```@contents
Pages = ["index.md","api.md","contrib.md"]
```

<!-- URLS -->
[infill-img]: docs/src/assets/infill_radius_white.gif
[infill-url]: https://faun.rc.fas.harvard.edu/saydjari/CloudCovErr/thr_test.mp4