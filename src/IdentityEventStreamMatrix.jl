
#TODO: Should also have mintime? Instead of assuming a zero start.
struct IdentityEventStreamVector{T <: Real} <: AbstractEventStreamVector{T}
    events :: Vector{T}
    δ :: T
    maxtime :: T
    nbins :: Int
    IdentityEventStreamVector(events, δ, maxtime) = new{typeof(δ)}(events, δ, maxtime, ceil(Int, maxtime/δ))
end


struct IdentityEventStreamMatrix{T <: Real, L} <: AbstractEventStreamMatrix{T, L}
    events :: Vector{Tuple{T, L}}
    labels :: Vector{L}
    δ :: T
    maxtime :: T
    nbins :: Int
    ncols :: Int
    IdentityEventStreamMatrix(eventstream, labels, δ, maxtime) = new{typeof(δ), eltype(labels)}(
        events,
        labels,
        δ, 
        maxtime, 
        ceil(Int, maxtime/δ),
        length(labels)
    )
end

# Identify set of labels from those occuring in eventstream
IdentityEventStreamMatrix(eventstream, δ :: Real, maxtime :: Real) = IdentityEventStreamMatrix(eventstream, unique(last(e) for e in eventstream), δ, maxtime)

# Identify maximum time from last occuring event.
IdentityEventStreamMatrix(eventstream, δ :: Real) = IdentityEventStreamMatrix(eventstream, δ, maximum(first(e) for e in eventstream))

Array(E :: IdentityEventStreamVector) = Vector(E)
function Vector(E :: IdentityEventStreamVector)
    v = zeros(Int, E.nbins)
    for e in events(E)
        v[whichbin(e, fineness(E))] += 1
    end
    return v
end

Array(E :: IdentityEventStreamMatrix) = Matrix(E)
function Matrix(E :: IdentityEventStreamMatrix) 
    M = zeros(Int, (E.nbins, E.ncols))
    label_order = Dict(l => i for (i, l) in enumerate(E.labels))
    for (t, l) in events(E)
        M[whichbin(t, E.δ), label_order[l]] += 1
    end
    return M
end


size(E :: AbstractEventStreamVector) = (ceil(Int, duration(E)/fineness(E)),)
size(E :: IdentityEventStreamMatrix) = (E.nbins, E.ncols)

function getindex(E :: IdentityEventStreamMatrix, I...)
    # get range of events which are in the I[1]th (time) bin
    binrange = searchsorted(events(E), binend(I[1], fineness(E)), by=(x) -> whichbin(x[1], fineness(E)))
    # of those, only count those which match the I[2]th label
    out = count(e-> e[2]==labels(E)[I[2]], events(E)[binrange])
    return out
end


function getindex(E :: AbstractEventStreamVector, i :: Int) 
    binrange = searchsorted(events(E), binend(i, fineness(E)), by=(x) -> whichbin(x, fineness(E)))
    return length(binrange)
end

function setindex(E :: IdentityEventStreamMatrix, I...)
    @error "Index alteration not supported"
end

function setindex!(E :: AbstractEventStreamVector, v, i :: Int)
    @error "Index alteration not supported"
end


function findall(E :: IdentityEventStreamVector)
    return [whichbin(e, fineness(E)) for e in events(E)]
end

function zeros(::Type{IdentityEventStreamVector}, δ, maxtime)
    return IdentityEventStreamVector(Float64[], δ, maxtime)
end


duration(E :: AbstractEventStreamVector) = E.maxtime
duration(E :: AbstractEventStreamMatrix) = E.maxtime

fineness(E :: AbstractEventStreamMatrix) = E.δ
fineness(E :: AbstractEventStreamVector) = E.δ

events(E :: AbstractEventStreamMatrix) = E.events
events(E :: AbstractEventStreamVector) = E.events

labels(E :: AbstractEventStreamMatrix) = E.labels

