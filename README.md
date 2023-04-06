# NetworkHistogram

[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![codecov](https://codecov.io/gh/SDS-EPFL/NetworkHistogram.jl/branch/main/graph/badge.svg?token=CT0HIA66V1)](https://codecov.io/gh/SDS-EPFL/NetworkHistogram.jl)
[![CI](https://github.com/SDS-EPFL/NetworkHistogram.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SDS-EPFL/NetworkHistogram.jl/actions/workflows/CI.yml)


Implementation of the network histogram for graphon estimation from the paper [Network histograms and universality of blockmodel approximation](https://doi.org/10.1073/pnas.1400374111) by Sofia C. Olhede and Patrick J. Wolfe.


## Installation

```julia
Pkg.add("NetworkHistogram")
```

## Usage

We fit the estimator using [`graphhist`](@ref graphhist) and then extract the estimated graphon matrix and node labels.

```julia
using NetworkHistogram
A = rand(0:1, 100, 100)

# approximate the graphon with a network histogram
hist = graphhist(A)

# get the graphist structure
estimate = hist.graphist

# get the estimated graphon matrix
sbm_matrix = estimate.θ

# get the estimated node labels
node_labels = estimate.node_labels
```

You can control the optimization process by modifying the rules used in the optimization. Check out the docs for more information.