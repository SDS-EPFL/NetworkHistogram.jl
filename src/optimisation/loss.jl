# default score is the log likelihood
function score(a::Assignment, g::Observations)
    return log_likelihood(a, g) / binomial(number_nodes(a), 2)
end
