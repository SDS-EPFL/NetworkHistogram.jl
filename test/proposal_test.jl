@testset "Proposal" begin
    A = [0 1 1 1 0 0 1 0
         1 0 1 1 0 0 0 0
         1 1 0 0 0 0 0 0
         1 1 0 0 0 0 0 1
         0 0 0 0 0 1 1 1
         0 0 0 0 1 0 1 1
         1 0 0 0 1 1 0 0
         0 0 0 1 1 1 0 0]
    h = 0.5
    node_labels = [1, 1, 1, 1, 2, 2, 2, 2]
    group_size = NetworkHistogram.GroupSize(8, 4)
    swap = (2, 5)
    assignment = NetworkHistogram.Assignment(A, node_labels, group_size)
    proposal = deepcopy(assignment)
    NetworkHistogram.make_proposal!(proposal, assignment, swap, A)
    reference_proposal = NetworkHistogram.Assignment(A, [1, 2, 1, 1, 1, 2, 2, 2],
                                                     group_size)

    @testset "update labels" begin
        @test proposal.node_labels[swap[1]] == reference_proposal.node_labels[swap[1]] == 2
        @test proposal.node_labels[swap[2]] == reference_proposal.node_labels[swap[2]] == 1
    end

    @testset "update realized edges" begin
        @test proposal.realized[1, 2] == reference_proposal.realized[1, 2] == 8
        @test proposal.realized[2, 1] == reference_proposal.realized[2, 1] == 8
        @test proposal.realized[1, 1] == reference_proposal.realized[1, 1] == 2
        @test proposal.realized[2, 2] == reference_proposal.realized[2, 2] == 2
    end

    @testset "fast likelihood update" begin
        # inside each group likelihood contribution
        theoretical_after_update = 2 * (2 * log(2 / 6) + log(4 / 6) * 4)
        # between group likelihood contribution
        theoretical_after_update += 8 * log(8 / 16) * 2
        @test proposal.likelihood[1] == theoretical_after_update / 8 ==
              reference_proposal.likelihood[1]
    end
end
