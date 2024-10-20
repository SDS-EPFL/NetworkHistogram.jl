struct AggregData{C, R, E, D, F}
    counts::AbstractMatrix{C}
    realized::AbstractMatrix{R}
    estimated_theta::AbstractMatrix{E}
    A::D
    log_likelihood::F
end
