function bootstrap(
        statistic::Function, data::AbstractMatrix, model::BlockModel,
        sampling::BootstrapSampling)
    t0 = tx(statistic(data))
    m = nrun(sampling)
    t1 = zeros_tuple(t0, m)
    data1 = copy(data)
    for i in 1:m
        draw_and_fill!(data1, model)
        for (j, t) in enumerate(tx(statistic(data1)))
            t1[j][i] = t
        end
    end
    return ParametricBootstrapSample(t0, t1, statistic, data, model, sampling)
end
