module EventStreamMatrix

export binstart, binend, binmid
export nextbinmid, prevbinmid
export whichbin
export FirstOrderBSplineEventStreamMatrix
export IdentityEventStreamMatrix, IdentityEventStreamVector

import Base: findall, size, getindex, setindex!, length, zeros
import SparseArrays: findnz

using LinearAlgebra
using BSplines

abstract type AbstractEventStreamMatrix{T, S} <: AbstractMatrix{T} end
abstract type AbstractEventStreamVector{T} <: AbstractVector{T} end
greet() = "Hell"

include("utils.jl")
include("IdentityEventStreamMatrix.jl")
include("BSplineEventStream.jl")

end # module
