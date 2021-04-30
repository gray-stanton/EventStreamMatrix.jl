using EventStreamMatrix
using LinearAlgebra
using Test

@testset "All tests" begin
    
    @testset "IdentityEventStreamVector" begin
        @testset "Binning" begin
        @test whichbin(0.0, 1.0) == 1
        @test whichbin(0.5, 0.5) == 1
        @test whichbin(0.6, 0.5) == 2
        @test binend(1, 0.5) == 0.5
        @test binend(2, 0.5) == 1.0
        @test binstart(1, 0.5) == 0.0
        @test binstart(2, 0.5) == 0.5
        @test binmid(1, 0.5) == 0.25
        @test binmid(2, 0.5) == 0.75
        @test nextbinmid(0.75, 0.5) == 0.75
        @test nextbinmid(0.00, 0.5) == 0.25
        @test nextbinmid(0.85, 0.5) == 1.25
        @test nextbinmid(0.65, 0.5) == 0.75
        # Left-open, right-close except for first bin
        @test whichbin(binend(1, 0.5), 0.5) == 1
        @test whichbin(binend(2, 0.5), 0.5) == 2
        @test whichbin(binstart(1, 0.5), 0.5) == 1
        @test whichbin(binstart(2, 0.5), 0.5) == 1
        @test whichbin(binstart(3, 0.5), 0.5) == 2
        @test whichbin(nextbinmid(0.85, 0.5), 0.5) == 3
        @test whichbin(nextbinmid(0.75, 0.5), 0.5) == 2
        @test whichbin(nextbinmid(0.65, 0.5), 0.5) == 2
        end

        eventtimes = [2.0, 3.0, 3.25, 3.5, 4.5]
        maxtime = 5.0
        δ = 0.5
        E1 = IdentityEventStreamVector(eventtimes, δ, maxtime)
        @testset "Creation" begin
            @test E1.nbins == 10
            @test duration(E1) == 5.0
            @test fineness(E1) == 0.5
            @test events(E1) == eventtimes
            @test events(zeros(IdentityEventStreamVector, 1.0, 10.0)) == Float64[]
        end
        @testset "Conversion" begin
            binned_events = [0, 0, 0, 1, 0, 1, 2, 0, 1, 0]
            @test Vector(E1) == binned_events
            @test findall(E1) == [4, 6, 7, 7, 9]
        end
        @testset "AbstractVector Interface" begin
           @test size(E1) == (10,)
           @test E1[1] == 0
           @test E1[4] == 1
           @test E1[7] == 2
        end

    end

    @testset "IdentityEventStreamMatrix" begin
            eventtimes = [2.0, 3.0, 3.25, 3.5, 4.5]
            labelled = ["a", "b", "a", "b", "a"]
            labs = ["a", "b"]
            eventstream = collect(zip(eventtimes, labelled))
            maxtime = 5.0
            δ = 0.5
            E1 = IdentityEventStreamMatrix(eventstream, labs, δ, maxtime)
        @testset "Creation" begin
            E2 = IdentityEventStreamMatrix(eventstream, δ, maxtime)
            E3 = IdentityEventStreamMatrix(eventstream, δ)
            #@test labels(E1) == labels(E2)
            @test duration(E2) == 5.0
            @test size(E1)[2] == 2
        end

        @testset "Conversion" begin
            @test Matrix(E1) == [0 0; 0 0; 0 0; 1 0; 0 0; 0 1; 1 1; 0 0; 1 0; 0 0]
            @test sparse(E1) == false
        end
    end

    @testset "BSplineEventStreamMatrix" begin
        eventtimes1 = [0.1, 0.6, 1.0, 1.25]
        labs1 = ["a", "a", "a", "a"]
        labs2 = ["a", "b", "a", "b"]
        events1 = collect(zip(eventtimes1, labs1))
        events2 = collect(zip(eventtimes1, labs2))
        labels1 = ["a"]
        labels2 = ["a", "b"]
        δ = 0.1
        maxtime = 1.5
        splineorder = 2
        breakpoints = [0.0, 0.01, 0.05, 0.1, 0.2]
        E1 = FirstOrderBSplineEventStreamMatrix(events1, labels1, δ, maxtime, splineorder, breakpoints)
        E2 = FirstOrderBSplineEventStreamMatrix(events2, labels2, δ, maxtime, splineorder, breakpoints)
        @testset "Construction" begin
            @test E1.ncols == 5
            @test E2.ncols == 10
            @test length(E2.basis) == 5
            @test E1.memory == 0.2
        end
        @testset "Interface" begin
            @test size(E1) == (15, 5)
            @test getindex(E1, 15, 1) == 0
            @test getindex(E1, 2, 1) == 0.0
            @test getindex(E1, 2, 3) >= 0.9999
            @test getindex(E2, 7, 3) == 0.0
            @test (Matrix(E1)[2, 3] - getindex(E1, 2, 3)) <= 0.0001
        end

        @testset "Multiplication" begin
            y = collect(1.0:15.0)
            b = collect(1.0:5.0)
            W1 = collect(5.0:-1.0:1.0)
            W2 = collect(1.0:15.0)
            G = WeightedGramMatrix(E1, W2)
            eventtimes = sort(rand(Float64, 100))
            E3 = FirstOrderBSplineEventStreamMatrix(collect(zip(eventtimes, repeat(["a"], 100))),
                 ["a"], 0.01, 1.0, 4, [0.0, 0.01, 0.02, 0.03, 0.04, 0.05])
            E4 = FirstOrderBSplineEventStreamMatrix(collect(zip(eventtimes, repeat(["a", "b"], 50))),
                 ["a", "b"], 0.01, 1.0, 4, [0.0, 0.01, 0.02, 0.03, 0.04, 0.05])
            ws = zeros(5)
            W3 = rand(Float64, 100)
            @test abs(sum((Matrix(E1) * (W1 .* b)) - XWb(E1, W1, b))) <= 0.0001
            @test abs(sum(Matrix(E1)' * (W2 .* y) - XtWy(E1, W2, y))) <= 0.0001
            @test abs(sum(Matrix(E1)' * diagm(W2) * Matrix(E1))) <= 0.0001
            @test abs(sum(Matrix(E1)' * diagm(W2) * Matrix(E1) * b)) <= 0.0001
            @test sum(abs.(Matrix(E3)'* Matrix(E3) - XtWX(E3, ones(E3.nbins)))) <= 0.0001
            @test sum(abs.(Matrix(E4)'* Matrix(E4) - XtWX(E4, ones(E4.nbins)))) <= 0.0001
            @test sum(abs.(Matrix(E4)' * Matrix(E4))*collect(1:16.0) - XtWXb(E4, ones(E4.nbins), collect(1:16.0))) <= 0.0001
            @test sum(abs.(Matrix(E4)'*diagm(W3) * Matrix(E4) - XtWX(E4, W3))) <= 0.0001
            @test sum(abs.(Matrix(E4)'*diagm(W3) * Matrix(E4)*collect(1:16.0) - XtWXb(E4, W3, collect(1:16.0)))) <= 0.0001
        end
    end

    @testset "BSplineDiscEventStreamMatrix" begin
        eventtimes = 100.0 *rand(Float64, 100)
        labs = repeat(1:2, 50)
        evns = collect(zip(eventtimes, labs))
        sort!(evns)
        E = FirstOrderDiscBSplineEventStreamMatrix(evns, ["a", "b"], 0.1, 100.0, 4, [0.0, 2.0, 4.0], true)
    end
end

