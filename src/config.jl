struct NetworkHistConfig{T1, T2<:NodeSwapRule, T3, T4}
    starting_assignment::T1
    swapping_rule::T2
    accept_rule::T3
    stopping_rule::T4
end