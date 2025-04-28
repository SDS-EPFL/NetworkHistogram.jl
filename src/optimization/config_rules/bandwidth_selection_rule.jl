abstract type KSelectionRule end
struct OracleK <: KSelectionRule
    K::Int
end

struct OracleH <: KSelectionRule
    H::Int
end



"""
    select_number_node_per_block(g::Observations, rule::KSelectionRule)

How to select the number of blocks `K` for the BlockModel model.

# Implemented rules
- `OracleK(K::Int)`: Use the oracle number of blocks `K`.
- `OracleH(H::Int)`: Use the oracle number of nodes per block `H`.

!!! info
    - The number of blocks `K` should be at most `n/2` where `n` is the number of nodes in
        the graph.
"""
select_number_node_per_block

function select_number_node_per_block(g, rule::OracleH)
    if rule.H > number_nodes(g) ÷ 2
        throw(ArgumentError("The number of nodes per block $(rule.H) is too large for the \
        number of nodes $(number_nodes(g)), it should be at most $(number_nodes(g)÷2)"))
    end
    if rule.H <= 1
        throw(ArgumentError("The number of nodes per block $(rule.H) is too small, it should \
        be at least 2"))
    end
    return rule.H
end

function select_number_node_per_block(g, rule::OracleK)
    nodes_per_block = number_nodes(g) ÷ rule.K
    return select_number_node_per_block(g, OracleH(nodes_per_block))
end
