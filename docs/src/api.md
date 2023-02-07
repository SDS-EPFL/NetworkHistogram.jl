```@contents
Pages = ["api.md"]
Depth = 1
```

# Assignment 

```@autodocs
Modules = [NetworkHistogram]
Pages   = ["assignment.jl", "group_numbering.jl"]
```


# Proposal

```@autodocs
Modules = [NetworkHistogram]
Pages   = ["proposal.jl"]
```


# Starting node labels

!!! note
    The groups will be of size `floor(h * n)` where `n` is the number of nodes if `h` is a
    float. If `h` is an integer, the groups will be of size `h`. The last group may be
    smaller if `n` is not perfectely divisible by the group size.


```@autodocs
Modules = [NetworkHistogram]
Pages   = ["starting_assignment_rule.jl"]
```