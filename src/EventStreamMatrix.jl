module EventStreamMatrix

export binstart, binend, binmid
export nextbinmid, prevbinmid
export whichbin
export FirstOrderBSplineEventStreamMatrix
export IdentityEventStreamMatrix, IdentityEventStreamVector
export WeightedGramMatrix
export duration, fineness, events
export zeros, size, getindex, setindex!, length, findall
export Matrix, Array, Vector
export XtWX, XtWXb, XtWy, XWb
export XWb!, XtWy!, XtWX!, XtWXb!
export mul!
export WeightedNormGramMatrix
export FirstOrderDiscBSplineEventStreamMatrix
export memorylengths_away

export AbstractEventStreamMatrix, AbstractEventStreamVector

import Base: findall, size, getindex, setindex!, length, zeros, Matrix, Array, Vector
import SparseArrays: findnz
import LinearAlgebra: mul!

using LinearAlgebra
using BSplines
using DataStructures

abstract type AbstractEventStreamMatrix{T, S} <: AbstractMatrix{T} end
abstract type AbstractEventStreamVector{T} <: AbstractVector{T} end
greet() = "Hellp"

include("utils.jl")
include("IdentityEventStreamMatrix.jl")
include("BSplineEventStream.jl")

end # module
