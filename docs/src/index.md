# NetworkHistogram.jl

Implementation of the network histogram for graphon estimation from the paper [Network histograms and universality of blockmodel approximation](https://doi.org/10.1073/pnas.1400374111) by Sofia C. Olhede and Patrick J. Wolfe.


## Installation

```julia
Pkg.add("NetworkHistogram")
```

## Usage

We fit the estimator using [`graphhist`](@ref graphhist) and then extract the estimated graphon matrix and node labels.

```julia
using NetworkHistogram, LinearAlgebra

A = Symmetric(rand(0:1, 100, 100))
A[diagind(A)] .= 0

# approximate the graphon with a network histogram
hist = graphhist(A)

# get the graphist structure
estimate = hist.graphhist

# get the estimated graphon matrix
sbm_matrix = estimate.θ

# get the estimated node labels
node_labels = estimate.node_labels
```

You can control the optimization process by modifying the rules used in the optimization. Check out [Optimization hyper-parameters](@ref) for more information.