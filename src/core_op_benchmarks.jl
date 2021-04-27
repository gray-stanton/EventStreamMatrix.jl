using ArgParse
using EventStreamMatrix
using BenchmarkTools
using LinearAlgebra
using SparseArrays

import Random: seed!

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin # todo change all to be varargs '+'
        "-m"
            help="Number of events per source"
            default=100
            arg_type=Int
        "-n"
            help="Number of sources"
            default=1
            arg_type=Int
        "-t"
            help="Maximum time to simulate to"
            default=1.0
            arg_type=Float64
        "--fineness"
            help="Level of discretization"
            default=0.1
            arg_type=Float64
        "--seed"
            help="Random seed"
            default=234
            arg_type=Int
        "--seconds"
            help="Number of seconds to benchmark for"
            default=10.0
            arg_type=Float64
        "--breakpoints"
            help="Breakpoints for bspline"
            nargs='+'
            arg_type=Float64
            default=[0.0, 0.4, 0.8]
    end
    return parse_args(s)
end

function main()
    BLAS.set_num_threads(1)
    args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in args
        println("  $arg  =>  $val")
    end
    seed!(args["seed"])
    eventtimes = args["t"] * rand(Float64, args["m"] * args["n"])
    labels = ["s$i" for i in 1:args["n"]]
    labs = repeat(labels, floor(Int, length(eventtimes)/args["n"]))
    eventstream = collect(zip(eventtimes, labs))
    sort!(eventstream)
    E = FirstOrderBSplineEventStreamMatrix(eventstream, labels, args["fineness"], args["t"], 4, args["breakpoints"])
    b = randn(Float64, size(E)[2])
    y = randn(Float64, size(E)[1])
    W1 = rand(Float64, size(E)[1])
    W1d = Diagonal(W1)
    W2 = rand(Float64, size(E)[2])
    W2d = Diagonal(W2)
    dest1 = zeros(Float64, size(E)[1])
    dest2 = zeros(Float64, size(E)[2])
    dest3 = zeros(Float64, size(E)[2], size(E)[2])
    X = Matrix(E)
    S = sparse(X)
    #XWb trials
    print("XWb trials begin")
    xwb_dense = @benchmark $X * $W2d * $b seconds=args["seconds"]
    xwb_sparse = @benchmark $S * $W2d * $b seconds=args["seconds"]
    xwb_event = @benchmark XWb!($dest1, $E, $W2, $b) seconds=args["seconds"]

    #XtWy trials
    print("XtWy trials begin")
    xtwy_dense = @benchmark $X' * $W1d * $y seconds=args["seconds"]
    xtwy_sparse = @benchmark $S' * $W1d * $y seconds=args["seconds"]
    xtwy_event = @benchmark XtWy!($dest1, $E, $W1, $y) seconds=args["seconds"]

    #XtWX trials
    print("XtWX trials begin")
    xtwx_dense = @benchmark $X' * $W1d * $X seconds=args["seconds"]
    xtwx_sparse = @benchmark $S' * $W1d * $S seconds=args["seconds"]
    xtwx_event = @benchmark XtWX!($dest3, $E, $W1) seconds=args["seconds"]

    # XtWXb trials
    print("XtWXb trials begin")
    xtwxb_dense = @benchmark $X' * $W1d * $X * $b seconds=args["seconds"]
    xtwxb_sparse = @benchmark $S' * $W1d * $S * $b seconds=args["seconds"]
    xtwxb_event = @benchmark XtWXb!($dest2, $E, $W1, $b) seconds=args["seconds"]

    names = ["xwb_dense", "xwb_sparse", "xwb_event", "xtwy_dense", "xtwy_sparse", "xtwy_event", "xtwx_dense",
             "xtwx_spase", "xtwx_event", "xtwxb_dense", "xtwxb_sparse", "xtwxb_event"]
    trials=[xwb_dense, xwb_sparse, xwb_event, xtwy_dense, xtwy_sparse, xtwy_event, xtwx_dense, 
            xtwx_spase, xtwx_event, xtwxb_dense, xtwxb_sparse, xtwxb_event]
    for (n,t) in zip(trials, trials)
        println(n)
        println(t)
    end
end

main()
