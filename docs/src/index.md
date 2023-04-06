# NetworkHistogram.jl

Implementation of the network histogram for graphon estimation from the paper [Network histograms and universality of blockmodel approximation](https://doi.org/10.1073/pnas.1400374111) by Sofia C. Olhede and Patrick J. Wolfe.


## Installation

```julia
Pkg.add("NetworkHistogram")
```

## Usage

```julia
using NetworkHistogram
A = rand(0:1, 100, 100)

# approximate the graphon with a network histogram
hist = network_histogram(A)

# get the graphist structure
estimate = hist.graphist

# get the estimated graphon matrix
sbm_matrix = estimate.θ

# get the estimated node labels
node_labels = estimate.node_labels
```