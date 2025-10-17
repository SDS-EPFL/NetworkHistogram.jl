"""
Profiling helper for NetworkHistogram optimization.

This script helps identify performance bottlenecks using Julia's Profile module.

Usage:
    julia --project=. benchmark/profile_optimization.jl [scenario]

Scenarios:
    swap        - Profile single swap operations
    optimize    - Profile full optimization run
    components  - Profile individual components

The script generates profiling data and can display it as:
- Text output (default)
- FlameGraph (requires ProfileView.jl or PProf.jl)
"""

using Profile
using Random
using StatsBase
using StaticArrays
using NetworkHistogram

# Helper functions to create test networks
function create_test_sbm_bernoulli(n_groups::Int, n_nodes::Int; seed = 42)
    Random.seed!(seed)
    d = NetworkHistogram.Bernoulli(0.5)
    sbm = NetworkHistogram.BlockModel(n_groups, d)

    for g1 in 1:n_groups
        for g2 in g1:n_groups
            p = 0.1 + 0.7 * rand()
            sbm[g1, g2] = NetworkHistogram.Bernoulli(p)
        end
    end

    base_size = n_nodes ÷ n_groups
    remainder = n_nodes % n_groups
    sizes = fill(base_size, n_groups)
    sizes[1:remainder] .+= 1
    labels = StatsBase.inverse_rle(1:n_groups, sizes)
    A = NetworkHistogram.sample(sbm, labels)
    return A, labels, d
end

function create_test_sbm_categorical(
        n_groups::Int, n_nodes::Int, n_categories::Int; seed = 42)
    Random.seed!(seed)
    ps = SVector{n_categories}(fill(1 / n_categories, n_categories))
    d = NetworkHistogram.Cat(ps)
    sbm = NetworkHistogram.BlockModel(n_groups, d)

    for g1 in 1:n_groups
        for g2 in g1:n_groups
            probs = rand(n_categories)
            probs ./= sum(probs)
            sbm[g1, g2] = NetworkHistogram.Cat(SVector{n_categories}(probs))
        end
    end

    base_size = n_nodes ÷ n_groups
    remainder = n_nodes % n_groups
    sizes = fill(base_size, n_groups)
    sizes[1:remainder] .+= 1
    labels = StatsBase.inverse_rle(1:n_groups, sizes)
    A = NetworkHistogram.sample(sbm, labels)
    return A, labels, d
end

function profile_swap_operations(network_type = :bernoulli, n = 200, k = 3)
    println("Setting up $(network_type) network (n=$n, k=$k)...")

    if network_type == :bernoulli
        A, labels, d = create_test_sbm_bernoulli(k, n)
    else
        A, labels, d = create_test_sbm_categorical(k, n, 4)
    end

    edgelist = NetworkHistogram.EdgeList(A)
    assignment = NetworkHistogram.Assignment(labels, edgelist, NetworkHistogram.Dist(d))
    swap = NetworkHistogram.make_swap(assignment, (1, n))

    # Warm up
    println("Warming up...")
    for i in 1:100
        NetworkHistogram.apply_swap!(assignment, swap)
        NetworkHistogram.revert_swap!(assignment, swap)
    end

    # Profile
    println("Profiling swap operations (5000 iterations)...")
    Profile.clear()
    @profile begin
        for i in 1:5000
            NetworkHistogram.apply_swap!(assignment, swap)
            NetworkHistogram.revert_swap!(assignment, swap)
        end
    end

    return true
end

function profile_full_optimization(
        network_type = :bernoulli, n = 200, k = 3, max_iter = 10_000)
    println("Setting up $(network_type) network (n=$n, k=$k)...")

    if network_type == :bernoulli
        A, labels, d = create_test_sbm_bernoulli(k, n)
    else
        A, labels, d = create_test_sbm_categorical(k, n, 4)
    end

    initial_labels = rand(1:k, n)
    params = NetworkHistogram.GreedyParams(
        max_iter,
        NetworkHistogram.RandomNodeSwap(),
        NetworkHistogram.Strict(),
        NetworkHistogram.PreviousBestValue(max_iter),
        false
    )

    # Warm up
    println("Warming up...")
    test_params = NetworkHistogram.GreedyParams(
        100,
        NetworkHistogram.RandomNodeSwap(),
        NetworkHistogram.Strict(),
        NetworkHistogram.PreviousBestValue(50),
        false
    )
    NetworkHistogram.nethist(A, d, copy(initial_labels), test_params)

    # Profile
    println("Profiling full optimization ($max_iter iterations)...")
    Profile.clear()
    @profile NetworkHistogram.nethist(A, d, initial_labels, params)

    return true
end

function profile_components(n = 200, k = 3)
    println("Setting up network (n=$n, k=$k)...")
    A, labels, d = create_test_sbm_bernoulli(k, n)
    edgelist = NetworkHistogram.EdgeList(A)

    # Profile Assignment creation
    println("\nProfiling Assignment creation...")
    Profile.clear()
    @profile begin
        for i in 1:1000
            NetworkHistogram.Assignment(labels, edgelist, NetworkHistogram.Dist(d))
        end
    end
    println("Results for Assignment creation:")
    Profile.print(maxdepth = 15)

    # Profile EdgeList creation
    println("\nProfiling EdgeList creation...")
    Profile.clear()
    @profile begin
        for i in 1:1000
            NetworkHistogram.EdgeList(A)
        end
    end
    println("Results for EdgeList creation:")
    Profile.print(maxdepth = 15)

    # Profile log-likelihood computation
    assignment = NetworkHistogram.Assignment(labels, edgelist, NetworkHistogram.Dist(d))
    println("\nProfiling log-likelihood computation...")
    Profile.clear()
    @profile begin
        for i in 1:10000
            NetworkHistogram.loglikelihood(assignment)
        end
    end
    println("Results for log-likelihood computation:")
    Profile.print(maxdepth = 15)

    # Profile get_edges_in_groups
    println("\nProfiling get_edges_in_groups...")
    Profile.clear()
    @profile begin
        for i in 1:10000
            NetworkHistogram.get_edges_in_groups(assignment, 1, 2)
        end
    end
    println("Results for get_edges_in_groups:")
    Profile.print(maxdepth = 15)

    return false  # Don't print again at the end
end

function print_results()
    println("\n" * "="^70)
    println("Profile Results")
    println("="^70)
    println("\nTop functions by exclusive time:")
    Profile.print(maxdepth = 15)

    println("\n" * "="^70)
end

function print_help()
    println("""
    NetworkHistogram Profiling Helper
    ==================================

    Usage: julia --project=. dev/profile_optimization.jl [scenario] [options]

    Scenarios:
      swap             Profile swap operations (default)
      swap-bernoulli   Profile Bernoulli swap operations
      swap-categorical Profile Categorical swap operations
      optimize         Profile full optimization run
      components       Profile individual components
      help             Show this message

    Options:
      --n=N            Number of nodes (default: 200)
      --k=K            Number of groups (default: 3)
      --iter=N         Number of iterations for optimization (default: 10000)

    Examples:
      julia dev/profile_optimization.jl swap
      julia dev/profile_optimization.jl swap-categorical --n=500
      julia dev/profile_optimization.jl optimize --iter=5000
      julia dev/profile_optimization.jl components

    The profiling results will show:
      - Which functions consume the most time
      - Call counts for each function
      - Memory allocation patterns
      - Call stack visualization (with flamegraph viewer)

    Tips:
      - Focus on functions with high "exclusive" time
      - Look for unexpected allocations
      - Check for type instabilities
      - Use flamegraph for visual exploration
    """)
end

function parse_args(args)
    options = Dict(
        :n => 200,
        :k => 3,
        :iter => 10_000
    )

    for arg in args
        if startswith(arg, "--n=")
            options[:n] = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--k=")
            options[:k] = parse(Int, split(arg, "=")[2])
        elseif startswith(arg, "--iter=")
            options[:iter] = parse(Int, split(arg, "=")[2])
        end
    end

    return options
end

function main()
    if length(ARGS) == 0 || ARGS[1] in ["help", "-h", "--help"]
        print_help()
        return
    end

    scenario = ARGS[1]
    options = parse_args(ARGS[2:end])

    should_print = if scenario == "swap" || scenario == "swap-bernoulli"
        profile_swap_operations(:bernoulli, options[:n], options[:k])
    elseif scenario == "swap-categorical"
        profile_swap_operations(:categorical, options[:n], options[:k])
    elseif scenario == "optimize"
        profile_full_optimization(:bernoulli, options[:n], options[:k], options[:iter])
    elseif scenario == "components"
        profile_components(options[:n], options[:k])
    else
        println("Error: Unknown scenario '$scenario'")
        print_help()
        return
    end

    if should_print
        print_results()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
