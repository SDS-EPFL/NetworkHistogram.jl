"""Utils to read .gt files

Inspired from the library Erdos.jl from CarloLucibello (precisely the file
`src/persistence.jl`)
"""

const start_gt_format = "⛾ gt"

function readgt_adj!(io::IO, adj, ::Type{T}) where {T}
    n = size(adj, 1)
    for i in 1:n
        k = read(io, UInt64)
        for _ in 1:k
            j = read(io, T) + 1
            adj[i, j] = 1
            adj[j, i] = 1
        end
    end
end

function readgt(io::IO)
    @assert String(read(io, 6))==start_gt_format "gt file not correctly formatted"
    ver = read(io, UInt8)  ## version
    indian = read(io, Bool)
    @assert indian == false
    lencomment = read(io, UInt64)
    read(io, lencomment)
    isdir = read(io, Bool)
    n = read(io, UInt64)
    T = minutype(n)
    if isdir
        @warn "Directed graphs are not supported, automatically converting to undirected by
        dropping the direction"
    end
    adj = zeros(Int, n, n)

    readgt_adj!(io, adj, T)
    return adj
end
