module EventStreamMatrix

import Base: findall, size, getindex, setindex!, length, zeros
import SparseArrays: findnz

using LinearAlgebra
using BSplines

abstract type AbstractEventStreamMatrix{T, S} <: AbstractMatrix{T} end
abstract type AbstractEventStreamVector{T} <: AbstractVector{T} end
end # module
