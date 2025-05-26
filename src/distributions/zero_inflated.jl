
struct ZeroInflated{D, F}
    dist::D
    proba_zero::F
end

function ZeroInflated(dist)
    return ZeroInflated(dist, 0.0)
end
