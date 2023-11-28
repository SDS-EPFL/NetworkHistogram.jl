var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API Reference","title":"API Reference","text":"Pages = [\"api.md\"]\nDepth = 1","category":"page"},{"location":"api/#NetworkHistogram","page":"API Reference","title":"NetworkHistogram","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Modules = [NetworkHistogram]\nPages   = [\"histogram.jl\",\"optimize.jl\"]","category":"page"},{"location":"api/#NetworkHistogram.GraphHist","page":"API Reference","title":"NetworkHistogram.GraphHist","text":"Network Histogram approximation [1].\n\nContains the estimated network histogram and the node labels.\n\nFields\n\nθ::Matrix{T}: Estimated stochastic block model parameters.\nnode_labels::Vector{Int}: Node labels for each node in the adjacency matrix used   to estimate the network histogram.\n\nReferences\n\n[1] - Olhede, Sofia C., and Patrick J. Wolfe. \"Network histograms and universality of blockmodel approximation.\" Proceedings of the National Academy of Sciences 111.41 (2014): 14722-14727.\n\n\n\n\n\n","category":"type"},{"location":"api/#NetworkHistogram._graphhist","page":"API Reference","title":"NetworkHistogram._graphhist","text":"_graphhist(A, record_trace=Val{true}(); h, maxitr, swap_rule, starting_assignment_rule, accept_rule, stop_rule)\n\nInternal version of graphhist which is type stable.\n\n\n\n\n\n","category":"function"},{"location":"api/#NetworkHistogram.graphhist-Tuple{Any}","page":"API Reference","title":"NetworkHistogram.graphhist","text":"graphhist(A; h = select_bandwidth(A), maxitr = 1000, swap_rule = RandomNodeSwap(),\nstarting_assignment_rule = RandomStart(), accept_rule = Strict(),\nstop_rule = PreviousBestValue(3), record_trace=true)\n\nComputes the graph histogram approximation.\n\nArguments\n\nA: adjacency matrix of a simple graph\nh: bandwidth of the graph histogram (number of nodes in a group or percentage (in [0,1]) of   nodes in a group)\nrecord_trace (optional): whether to record the trace of the optimization process and return   it as part of the output. Default is true.\n\nReturns\n\nnamed tuple with the following fields:\n\ngraphhist: the graph histogram approximation\ntrace: the trace of the optimization process (if record_trace is true)\nlikelihood: the loglikelihood of the graph histogram approximation\n\nExamples\n\njulia> A = [0 0 1 0 1 0 1 1 0 1\n             0 0 1 1 1 1 1 1 0 0\n             1 1 0 1 0 0 0 0 1 0\n             0 1 1 0 1 0 1 0 0 0\n             1 1 0 1 0 0 1 0 0 1\n             0 1 0 0 0 0 0 1 0 0\n             1 1 0 1 1 0 0 1 0 1\n             1 1 0 0 0 1 1 0 0 1\n             0 0 1 0 0 0 0 0 0 1\n             1 0 0 0 1 0 1 1 1 0]\njulia> out = graphhist(A);\njulia> graphist_approx = out.graphist\n...\njulia> trace = out.trace\nNetworkHistogram.TraceHistory{...}\n  :best_likelihood => 1 elements {Int64,Float64}\n  :proposal_likelihood => 5 elements {Int64,Float64}\n  :current_likelihood => 5 elements {Int64,Float64})\njulia> loglikelihood = out.likelihood\n-22.337057781338277\n\n\n\n\n\n","category":"method"},{"location":"api/#NetworkHistogram.graphhist_format_output-Tuple{Any, NetworkHistogram.TraceHistory}","page":"API Reference","title":"NetworkHistogram.graphhist_format_output","text":"graphhist_format_output(best, history)\n\nFormates the graphhist output depending on the type of history requested by the user.\n\n\n\n\n\n","category":"method"},{"location":"api/#NetworkHistogram.initialize-NTuple{4, Any}","page":"API Reference","title":"NetworkHistogram.initialize","text":"initialize(A, h, starting_assignment_rule, record_trace)\n\nInitialize the memory required for finding optimal graph histogram.\n\n\n\n\n\n","category":"method"},{"location":"api/#NetworkHistogram.oracle_bandwidth","page":"API Reference","title":"NetworkHistogram.oracle_bandwidth","text":"oracle_bandwidth(A, type = \"degs\", alpha = 1, c = min(4, sqrt(size(A, 1)) / 8))\n\nOracle bandwidth selection for graph histogram, using\n\nwidehath^*=left(2left(left(d^T dright)^+right)^2 d^T A d cdot hatm hatbright)^-frac12 hatrho_n^frac14\n\nwhere d is the vector of degree sorted in increasing order,hatrho_n is the empirical edge density, and  m, b are the slope and intercept fitted on dn2-csqrtnn2+csqrtn for some c.\n\n\n\n\n\n","category":"function"},{"location":"api/#Assignment","page":"API Reference","title":"Assignment","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Modules = [NetworkHistogram]\nPages   = [\"assignment.jl\", \"group_numbering.jl\"]","category":"page"},{"location":"api/#NetworkHistogram.compute_log_likelihood-NTuple{4, Any}","page":"API Reference","title":"NetworkHistogram.compute_log_likelihood","text":"compute_log_likelihood(number_groups, estimated_theta, counts, number_nodes)\n\nCompute the scaled log-likelihood in terms of communities:\n\nl(zA) = frac1n sum_g_1 = 1^G sum_g_2 geq g_1^G  left theta_g_1g_2 log(theta_g_1g_2) + (1 - theta_g_1g_2) log(1 - theta_g_1g_2) right cdot c_g_1g_2\n\nwhere c_g_1g_2 (theta_g_1g_2) is the number of possible edges (estimated probability) between communities g_1 and g_2, n is the number of nodes, and z_i  1 dots G is the community assignment of node i.\n\n\n\n\n\n","category":"method"},{"location":"api/#NetworkHistogram.compute_log_likelihood-Tuple{NetworkHistogram.Assignment}","page":"API Reference","title":"NetworkHistogram.compute_log_likelihood","text":"compute_log_likelihood(assignment::Assignment)\n\nCompute the scaled log-likelihood of the assignment.\n\n\tl(zA) = frac1nsumlimits_i=1^n sumlimits_ji^n  left A_ij log(hattheta_z_i z_j) + (1 - A_ij) log(1 - hattheta_z_i z_j) right\n\nwhere hattheta_ab is the estimated probability of an edge between communities a and b\n\n    hattheta_ab = fracsumlimits_ij A_ij mathbb1(z_i = a z_j = b) sumlimits_ij mathbb1(z_i = a z_j = b)\n\n\n\n\n\n","category":"method"},{"location":"api/#NetworkHistogram.GroupSize","page":"API Reference","title":"NetworkHistogram.GroupSize","text":"Array-like storage for the number of nodes in each group.\n\n\n\n\n\n","category":"type"},{"location":"api/#Proposal","page":"API Reference","title":"Proposal","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Modules = [NetworkHistogram]\nPages   = [\"proposal.jl\"]","category":"page"},{"location":"api/#NetworkHistogram.create_proposal!-Tuple{NetworkHistogram.GraphOptimizationHistory, Int64, NetworkHistogram.Assignment, NetworkHistogram.Assignment, Any, Any}","page":"API Reference","title":"NetworkHistogram.create_proposal!","text":"create_proposal!(history::GraphOptimizationHistory, iteration::Int, proposal::Assignment,\n                      current::Assignment, A, swap_rule)\n\nCreate a new proposal by swapping the labels of two nodes. The new assignment is stored in proposal. The swap is selected using the swap_rule function. The likelihood of the new proposal is stored in the history.\n\nwarning: Warning\nThe proposal assignment is modified in place to avoid unnecessary memory allocation.\n\n\n\n\n\n","category":"method"},{"location":"api/#NetworkHistogram.make_proposal!-Tuple{NetworkHistogram.Assignment, NetworkHistogram.Assignment, Tuple{Int64, Int64}, Any}","page":"API Reference","title":"NetworkHistogram.make_proposal!","text":"make_proposal!(proposal::Assignment, current::Assignment, swap::Tuple{Int, Int}, A)\n\nFrom the current assignment, create a new assignment by swapping the labels of the nodes specified in swap. The new assignment is stored in proposal.\n\n\n\n\n\n","category":"method"},{"location":"api/#NetworkHistogram.updateLL!-Tuple{NetworkHistogram.Assignment}","page":"API Reference","title":"NetworkHistogram.updateLL!","text":"updateLL!(proposal::Assignment)\n\nUpdate the likelihood of the proposal assignment based on its observed and estimated attributes.\n\n\n\n\n\n","category":"method"},{"location":"api/#NetworkHistogram.update_labels!-Tuple{NetworkHistogram.Assignment, Tuple{Int64, Int64}, NetworkHistogram.Assignment}","page":"API Reference","title":"NetworkHistogram.update_labels!","text":"update_labels!(proposal::Assignment, swap::Tuple{Int, Int}, current::Assignment)\n\nUpdate the labels of the nodes specified in swap in the proposal assignment.\n\n\n\n\n\n","category":"method"},{"location":"api/#NetworkHistogram.update_observed!-Tuple{NetworkHistogram.Assignment, Tuple{Int64, Int64}, Any}","page":"API Reference","title":"NetworkHistogram.update_observed!","text":"update_observed!(proposal::Assignment, swap::Tuple{Int, Int}, A)\n\nUpdate the observed and estimated attributes of the proposal assignment based on the swap of the nodes specified in swap.\n\nNOTE labels of the nodes before the swap\n\n\n\n\n\n","category":"method"},{"location":"#NetworkHistogram.jl","page":"Home","title":"NetworkHistogram.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Implementation of the network histogram for graphon estimation from the paper Network histograms and universality of blockmodel approximation by Sofia C. Olhede and Patrick J. Wolfe.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pkg.add(\"NetworkHistogram\")","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We fit the estimator using graphhist and then extract the estimated graphon matrix and node labels.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using NetworkHistogram, LinearAlgebra\n\nA = Symmetric(rand(0:1, 100, 100))\nA[diagind(A)] .= 0\n\n# approximate the graphon with a network histogram\nhist = graphhist(A)\n\n# get the graphist structure\nestimate = hist.graphhist\n\n# get the estimated graphon matrix\nsbm_matrix = estimate.θ\n\n# get the estimated node labels\nnode_labels = estimate.node_labels","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can control the optimization process by modifying the rules used in the optimization. Check out Optimization hyper-parameters for more information.","category":"page"},{"location":"internals/#Notes-on-optimisation","page":"Development","title":"Notes on optimisation","text":"","category":"section"},{"location":"internals/","page":"Development","title":"Development","text":"We have three different assignment variables:\ncurrent\nbest\nproposal\nWe need best because we might accept a proposal which is worse than the curent value.\nThis is dealt within accept_reject_update!().","category":"page"},{"location":"rules/#Optimization-hyper-parameters","page":"Optimization hyperparameters","title":"Optimization hyper-parameters","text":"","category":"section"},{"location":"rules/","page":"Optimization hyperparameters","title":"Optimization hyperparameters","text":"Here we discuss the different parameters that can be used to control the optimization process. The optimization greedily tries to find a good partition of the network. You can control the optimization process by setting the following parameters:","category":"page"},{"location":"rules/#Starting-node-labels","page":"Optimization hyperparameters","title":"Starting node labels","text":"","category":"section"},{"location":"rules/","page":"Optimization hyperparameters","title":"Optimization hyperparameters","text":"Modules = [NetworkHistogram]\nPages   = [\"starting_assignment_rule.jl\"]","category":"page"},{"location":"rules/#NetworkHistogram.initialize_node_labels","page":"Optimization hyperparameters","title":"NetworkHistogram.initialize_node_labels","text":"initialize_node_labels(A, h, starting_assignment_rule::StartingAssignment)\n\ninitialize node labels based on the starting_assignment_rule, and return a vector of node labels and a GroupSize object.\n\nImplemenented rules\n\nOrderedStart(): Sequentially assign nodes to groups based on the ordering of A.\nRandomStart(): Randomly assign nodes to groups.\nEigenStart(): Assign nodes to groups based on the second eigenvector of the normalized Laplacian.\n\n\n\n\n\n","category":"function"},{"location":"rules/","page":"Optimization hyperparameters","title":"Optimization hyperparameters","text":"note: Note\nThe groups will be of size floor(h * n) where n is the number of nodes if h is a float. If h is an integer, the groups will be of size h. The last group may be smaller if n is not exactly divisible by the group size.","category":"page"},{"location":"rules/#Swapping-rule","page":"Optimization hyperparameters","title":"Swapping rule","text":"","category":"section"},{"location":"rules/","page":"Optimization hyperparameters","title":"Optimization hyperparameters","text":"Modules = [NetworkHistogram]\nPages   = [\"swap_rule.jl\"]","category":"page"},{"location":"rules/#NetworkHistogram.select_swap","page":"Optimization hyperparameters","title":"NetworkHistogram.select_swap","text":"select_swap(node_assignment::Assignment, A, ::NodeSwapRule)\n\nSelects two nodes to swap based on the NodeSwapRule, the adjacency matrix A and the current assignment node_assignment.\n\nImplemented rules\n\nRandomNodeSwap(): Select two nodes at random.\n\n\n\n\n\n","category":"function"},{"location":"rules/#Acceptance-rule","page":"Optimization hyperparameters","title":"Acceptance rule","text":"","category":"section"},{"location":"rules/","page":"Optimization hyperparameters","title":"Optimization hyperparameters","text":"Modules = [NetworkHistogram]\nPages   = [\"accept_rule.jl\"]","category":"page"},{"location":"rules/#NetworkHistogram.accept_reject_update!","page":"Optimization hyperparameters","title":"NetworkHistogram.accept_reject_update!","text":"accept_reject_update!(history::GraphOptimizationHistory, iteration::Int,\n                      proposal::Assignment,\n                      current::Assignment, accept_rule::AcceptRule)\n\nReturn the updated current assignment based on the accept_rule.\n\nImplemented rules\n\nStrict(): Accept the proposal if it has a higher likelihood than the current assignment.\n\n\n\n\n\n","category":"function"},{"location":"rules/#Stopping-rule","page":"Optimization hyperparameters","title":"Stopping rule","text":"","category":"section"},{"location":"rules/","page":"Optimization hyperparameters","title":"Optimization hyperparameters","text":"Modules = [NetworkHistogram]\nPages   = [\"stop_rule.jl\"]","category":"page"},{"location":"rules/#NetworkHistogram.stopping_rule","page":"Optimization hyperparameters","title":"NetworkHistogram.stopping_rule","text":"stopping_rule(history, stop_rule::StopRule)\n\nReturns a Bool with true if we should stop the optimization based on the stop_rule.\n\nImplemented rules\n\nPreviousBestValue(k): Stop if the current iteration is k iterations away from the iteration with the best value.\n\n\n\n\n\n","category":"function"}]
}