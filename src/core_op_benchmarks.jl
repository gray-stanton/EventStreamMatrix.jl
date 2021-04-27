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
        "--breakpoints"
            help="Breakpoints for bspline"
            nargs="+"
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

end

main()
