# Optimization hyper-parameters

Here we discuss the different parameters that can be used to control the optimization process. The optimization greedily tries to find a good partition of the network. You can control the optimization process by setting the following parameters:

## Starting node labels

```@autodocs
Modules = [NetworkHistogram]
Pages   = ["starting_assignment_rule.jl"]
```

!!! note
    The groups will be of size `floor(h * n)` where `n` is the number of nodes if `h` is a
    float. If `h` is an integer, the groups will be of size `h`. The last group may be
    smaller if `n` is not exactly divisible by the group size.


## Swapping rule

```@autodocs
Modules = [NetworkHistogram]
Pages   = ["swap_rule.jl"]
```


## Acceptance rule

```@autodocs
Modules = [NetworkHistogram]
Pages   = ["accept_rule.jl"]
```

## Stopping rule

```@autodocs
Modules = [NetworkHistogram]
Pages   = ["stop_rule.jl"]
```