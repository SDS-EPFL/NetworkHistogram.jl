# Optimization hyper-parameters

Here we discuss the different parameters that can be used to control the optimization process. The optimization greedily tries to find a good partition of the network. You can control the optimization process by setting the following parameters:

## Starting node labels

```@docs; canonical=false
NetworkHistogram.initialize_node_labels
```

!!! note
    The groups will be of size `floor(h * n)` where `n` is the number of nodes if `h` is a
    float. If `h` is an integer, the groups will be of size `h`. The last group may be
    bigger if `n` is not exactly divisible by the group size.


## Swapping rule

```@docs; canonical=false
NetworkHistogram.select_swap
```


## Acceptance rule

```@docs; canonical=false
NetworkHistogram.accept_reject_update!
```

## Stopping rule

```@docs; canonical=false
NetworkHistogram.stopping_rule
```