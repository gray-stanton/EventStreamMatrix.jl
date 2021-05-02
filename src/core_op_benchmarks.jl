using ArgParse
using EventStreamMatrix
using BenchmarkTools
using LinearAlgebra
using SparseArrays
using DataFrames
using CSV

import Random: seed!

function f(x, w, b) 
    return x*w*b
end
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin # todo change all to be varargs '+'
        "-m"
            help="Number of events per source"
            nargs='+'
            default=[100]
            arg_type=Int
        "--mdensity"
            nargs='+'
            help="Density of events per source"
            arg_type=Float64
            default=Float64[]
        "-n"
            help="Number of sources"
            nargs='+'
            default=[1]
            arg_type=Int
        "-t"
            help="Maximum time to simulate to"
            default=[1.0]
            nargs='+'
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
        "--samples"
            help="Number of benchmark samples"
            default = 10
            arg_type=Int
        "--outfile"
            help="File to write results to"
            arg_type=String
            default="out.csv"
        "--breakpoints"
            help="Breakpoints for bspline"
            nargs='+'
            arg_type=Float64
            default=[0.0, 0.4, 0.8]
        "--skipdense"
            help="Whether or not to skip all dense multiplies"
            action=:store_true
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
    times = Float64[]
    mems = Float64[]
    ts = Float64[]
    ns = Int[]
    ms = Float64[]
    sizes = Float64[]
    nrows = Int[]
    ncols = Int[]
    sparsity = Float64[]
    trialnames = String[]
    for t in args["t"] 
        for n in args["n"]
            if length(args["mdensity"]) > 0
                miter = args["mdensity"]
            else 
                miter = args["m"]
            end
            for m in miter
                println("Begin t:$t n:$n m:$m")
                if length(args["mdensity"]) > 0
                    eventtimes = t * rand(Float64, ceil(Int, (m*t)*n))
                else
                    eventtimes = t * rand(Float64, m * n)
                end
                labels = ["s$i" for i in 1:n]
                labs = repeat(1:n, floor(Int, length(eventtimes)/n))
                eventstream = collect(zip(eventtimes, labs))
                sort!(eventstream)
                E = FirstOrderBSplineEventStreamMatrix(eventstream, labels, args["fineness"], t, 4, args["breakpoints"], true)
                b = randn(Float64, size(E)[2])
                y = randn(Float64, size(E)[1])
                W1 = rand(Float64, size(E)[1])
                W1d = Diagonal(W1)
                W2 = rand(Float64, size(E)[2])
                W2d = Diagonal(W2)
                dest1 = zeros(Float64, size(E)[1])
                dest2 = zeros(Float64, size(E)[2])
                dest3 = zeros(Float64, size(E)[2], size(E)[2])
                if !args["skipdense"]
                    X = Matrix(E)
                    S = sparse(X)
                end
                #XWb trials
                println("XWb trials begin")
                if !args["skipdense"]
                    xwb_dense = @benchmark $X * $W2d * $b seconds=args["seconds"]
                    xwb_sparse = @benchmark $S * $W2d * $b seconds=args["seconds"]
                end
                xwb_event = @benchmark XWb!($dest1, $E, $W2, $b) seconds=args["seconds"]

                #XtWy trials
                println("XtWy trials begin")
                if !args["skipdense"]
                    xtwy_dense = @benchmark $X' * $W1d * $y seconds=args["seconds"]
                    xtwy_sparse = @benchmark $S' * $W1d * $y seconds=args["seconds"]
                end
                xtwy_event = @benchmark XtWy!($dest2, $E, $W1, $y) seconds=args["seconds"]

                #XtWX trials
                println("XtWX trials begin")
                if !args["skipdense"]
                    xtwx_dense = @benchmark $X' * $W1d * $X seconds=args["seconds"]
                    xtwx_sparse = @benchmark $S' * $W1d * $S seconds=args["seconds"]
                end
                xtwx_event = @benchmark XtWX!($dest3, $E, $W1) seconds=args["seconds"]

                # XtWXb trials
                println("XtWXb trials begin")
                if !args["skipdense"]
                    xtwxb_dense = @benchmark $X' * $W1d * $X * $b seconds=args["seconds"]
                    xtwxb_sparse = @benchmark $S' * $W1d * $S * $b seconds=args["seconds"]
                end
                xtwxb_event = @benchmark XtWXb!($dest2, $E, $W1, $b) seconds=args["seconds"]
                if !args["skipdense"]
                    names = ["xwb_dense", "xwb_sparse", "xwb_event", "xtwy_dense", "xtwy_sparse", "xtwy_event", "xtwx_dense",
                            "xtwx_sparse", "xtwx_event", "xtwxb_dense", "xtwxb_sparse", "xtwxb_event"]
                    trials=[xwb_dense, xwb_sparse, xwb_event, xtwy_dense, xtwy_sparse, xtwy_event, xtwx_dense, 
                            xtwx_sparse, xtwx_event, xtwxb_dense, xtwxb_sparse, xtwxb_event]
                    memsizes = repeat([Base.summarysize(X)/2^20, Base.summarysize(S)/2^20, Base.summarysize(E)/2^20], 4)
                else
                    names = ["xwb_event", "xtwy_event", "xtwx_event", "xtwxb_event"]
                    trials=[xwb_event, xtwy_event, xtwx_event, xtwxb_event]
                    memsizes = repeat([Bas.summarysize(E)/2^20], 4)
                end
                trialnames = vcat(trialnames, names)
                sizes = vcat(sizes, memsizes)
                for (na,tr) in zip(names, trials)
                    mst = median(tr).time/1e6 #ms
                    memt = median(tr).memory/2^20 #MiB, total allocated
                    println("$t $n $m   $na: $mst ms")
                    push!(times, mst)
                    push!(mems, memt)
                    push!(ts, t)
                    push!(ms, m)
                    push!(ns, n)
                    push!(nrows, size(E)[1])
                    push!(ncols, size(E)[2])
                    if !args["skipdense"]
                        push!(sparsity, nnz(S)/length(S))
                    end
                end
            end
        end
    end
    print("Writing to $(args["outfile"])")
    df = DataFrame(:t => ts, :n => ns, :m => ms, :runtime => times, :memalloc => mems, :memsize => sizes,
        :trial => trialnames, :nrow => nrows, :ncol => ncols, :fillratio => sparsity
        )
    CSV.write(abspath(args["outfile"]), df)
end

main()
