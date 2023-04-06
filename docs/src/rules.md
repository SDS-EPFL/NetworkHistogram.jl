# Optimization parameters

Here we discuss the different parameters that can be used to control the optimization process. The optimization greedily tries to find a good partition of the network. You can control the optimization process by setting the following parameters:

## Starting node labels

```@autodocs
Modules = [NetworkHistogram]
Pages   = ["starting_assignment_rule.jl"]
```

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