struct FirstOrderBSplineEventStreamMatrix{T <: Real, L} <: AbstractEventStreamMatrix{T, L}
    events :: Vector{Tuple{T, L}}
    labels :: Vector{L}
    δ :: T
    maxtime :: T
    nbins :: Int
    ncols :: Int
    basis :: BSplineBasis{Vector{T}}
    memory :: T
    function FirstOrderBSplineEventStreamMatrix(
        events,
        labels,
        δ,
        maxtime,
        splineorder,
        breakpoints
    )
        newbasis = BSplineBasis(splineorder, breakpoints)
        return new{typeof(δ), eltype(labels)}(
           events,
           labels,
           δ, 
           maxtime, 
           ceil(Int, maxtime/δ),
           length(labels) * length(newbasis),
           newbasis,
           maximum(knots(newbasis))
    )
    end
end

size(E :: FirstOrderBSplineEventStreamMatrix) = (E.nbins, E.ncols)




function getindex(E :: FirstOrderBSplineEventStreamMatrix, I...)
    # get range of events which can influence it
    rellabel = E.labels[ceil(Int, I[2]/length(E.basis))]
    reltime = binmid(I[1], fineness(E))
    relspline_num = (I[2] -1) % length(E.basis) + 1
    relspline = E.basis[relspline_num]
    # only spikes in range [reltime - memory, reltime] should be included
    function memorylengths_away(t, reltime, memory)
        t == reltime - memory && return 0 # check exact boundary
        return ceil(Int, (nextfloat(t)-reltime)/memory)
    end
    binrange = searchsorted(events(E), prevfloat(reltime), by=(tup) -> memorylengths_away(tup[1], reltime, E.memory))
    tdiffs = [reltime - e[1] for e in events(E)[binrange] if e[2] == rellabel]
    out = sum(relspline.(tdiffs))
    return out
end

function setindex!(E :: FirstOrderBSplineEventStreamMatrix, I...)
    @error "Not implemented"
end

Array(E :: FirstOrderBSplineEventStreamMatrix) = Matrix(E)
function Matrix(E :: FirstOrderBSplineEventStreamMatrix) 
    #TODO make dense type adjustable
    out = zeros(Float64, E.nbins, E.ncols)
    #TODO avoid reallocating
    #dest = zeros(Float64, order(E.basis))
    label_order = Dict(l => i for (i, l) in enumerate(E.labels))
    starts = 1:nsplines:(E.ncols)
    for (t, l) in events(E)
        i = label_order[l]
        firstpoint = nextbinmid(t, E.δ)
        lastpoint = min(prevbinmid(t + E.memory, E.δ), prevbinmid(E.maxtime, E.δ))
        points = (firstpoint:δ:lastpoint) .- t
        if length(points) == 0
            # occurs if event in final half of last bin.
            continue
        end
        update_mat = basismatrix(E.basis, points)
        update_bins = whichbin(firstpoint, E.δ):1:whichbin(lastpoint, E.δ)
        update_cols = starts[i]:1:(starts[i]+length(E.basis)-1)
        if(length(update_bins) != length(points))
            @warn "Possible floating point error, due to events occuring exactly on bin boundary"
            len = min(length(update_bins), length(update_cols))
            update_bins=update_bins[1:len]
            update_mat = [1:len, :] 
        end
        out[update_bins, update_cols] += update_mat
    end
    return out
end


function XWb(E :: FirstOrderBSplineEventStreamMatrix{T}, W:: Vector{T}, b :: Vector{T}) where T
    out = zeros(eltype(b), E.nbins)
    nsplines = length(E.basis)
    fo_splines = [Spline(E.basis, W[k:k+nsplines-1].*b[k:(k+nsplines-1)])  for k in 1:nsplines:(E.ncols)]
    label_order = Dict(l => i for (i, l) in enumerate(E.labels))
    ws = zeros(eltype(b), order(E.basis))
    for (t, l) in events(E)
        i = label_order[l]
        firstpoint = nextbinmid(t, E.δ)
        lastpoint = min(prevbinmid(t + E.memory, E.δ), prevbinmid(E.maxtime, E.δ))
        points = (firstpoint:δ:lastpoint) .- t
        # Broadcast over points
        update_vec = fo_splines[i].(points; workspace=ws)
        update_bins = whichbin(firstpoint, E.δ):1:whichbin(lastpoint, E.δ)
        out[update_bins] += update_vec
    end
    return out
end



function XtWy(E :: FirstOrderBSplineEventStreamMatrix{T}, W, y :: Vector{T}) where T
    out = zeros(eltype(y), E.ncols)
    nsplines = length(E.basis)
    label_order = Dict(l => i for (i, l) in enumerate(E.labels))
    starts = 1:nsplines:(E.ncols)
    for (t, l) in events(E)
        i = label_order[l]
        firstpoint = nextbinmid(t, E.δ)
        lastpoint = min(prevbinmid(t + E.memory, E.δ), prevbinmid(E.maxtime, E.δ))
        points = (firstpoint:δ:lastpoint) .- t
        bmat = basismatrix(E.basis, points)
        update_bins = whichbin(firstpoint, E.δ):1:whichbin(lastpoint, E.δ)
        relW = W[update_bins]
        rely = y[update_bins]
        update_vec = transpose(bmat) * (relW .* rely)
        out[starts[i]:(starts[i] +nsplines -1)] += update_vec
    end
    return out
end


function XtWX(E:: FirstOrderBSplineEventStreamMatrix, W=ones(E.nbins))
    out = zeros(Float64, E.ncols, E.ncols)
    label_order = Dict(l => i for (i, l) in enumerate(E.labels))
    nsplines = length(E.basis)
    starts = 1:nsplines:(E.ncols)
    for (j, (t1, l1)) in enumerate(events(E))
        for (t2, l2) in events(E)[j:end]
            if t2 > t1 + E.memory
                break
            else
                i1 = label_order[l1]
                i2 = label_order[l2]
                firstpoint = max(nextbinmid(t1, E.δ), nextbinmid(t2, E.δ))
                lastpoint = min(prevbinmid(t1 + E.memory, E.δ),  prevbinmid(E.maxtime, E.δ))
                points_for_t1 = (firstpoint:δ:lastpoint) .- t1
                points_for_t2 = (firstpoint:δ:lastpoint) .- t2
                bmat_for_t1 = basismatrix(E.basis, points_for_t1)
                bmat_for_t2 = basismatrix(E.basis, points_for_t2)
                update_bins = whichbin(firstpoint, E.δ):1:whichbin(lastpoint, E.δ)
                relW = W[update_bins]
                # This is the "above" diagonal block, hence want the smaller of the two to transpose
                update_mat = i1 <= i2 ? bmat_for_t1' * diagm(relW) * bmat_for_t2 : bmat_for_t2' * diag(relW) * bmat_for_t1
                # Symmetric, update both blocks
                out[starts[i1]:(starts[i1]+nsplines-1), starts[i2]:(starts[i2] + nsplines-1)] += update_mat
                if i1 != i2
                    # On-diagonal blocks only need one update
                    out[starts[i2]:(starts[i2]+nsplines-1), starts[i1]:(starts[i1] + nsplines-1)] += transpose(update_mat)
                end
            end
        end
    end
    return out
end





function XtWXb(E :: FirstOrderBSplineEventStreamMatrix{T}, W, b :: Vector{T}) where T
    out = zeros(Float64, E.ncols)
    label_order = Dict(l => i for (i, l) in enumerate(E.labels))
    nsplines = length(E.basis)
    starts = 1:nsplines:(E.ncols)
    for (j, (t1, l1)) in enumerate(events(E))
        for (t2, l2) in events(E)[j:end]
            if t2 > t1 + E.memory
                break
            else
                i1 = label_order[l1]
                i2 = label_order[l2]
                firstpoint = max(nextbinmid(t1, E.δ), nextbinmid(t2, E.δ))
                lastpoint = min(prevbinmid(t1 + E.memory, E.δ),  prevbinmid(E.maxtime, E.δ))
                points_for_t1 = (firstpoint:δ:lastpoint) .- t1
                points_for_t2 = (firstpoint:δ:lastpoint) .- t2
                bmat_for_t1 = basismatrix(E.basis, points_for_t1)
                bmat_for_t2 = basismatrix(E.basis, points_for_t2)
                update_bins = whichbin(firstpoint, E.δ):1:whichbin(lastpoint, E.δ)
                relW = W[update_bins]
                update_mat = bmat_for_t1' * diagm(relW) * bmat_for_t2 
                update_vec1 = update_mat * b[starts[i1]:(starts[i1] + nsplines -1)]
                update_vec2 = transpose(update_mat) * b[starts[i2]:(starts[i2] + nsplines -1)]
                out[starts[i1]:(starts[i1] + nsplines -1)] += update_vec1
                if i1!= i2
                    # Avoid double-counting 
                    out[starts[i2]:(starts[i2] + nsplines -2)] += update_vec2
                end
            end
        end
    end
    return out
end





function convolve(A :: BSplineEventStreamMatrix{T, S}, β :: Vector{T},  δ :: T) where {T <: Real, S <: AbstractVector{T}}
    nbins = ceil(Int, A.maxtime/δ)
    dest = zeros(T, nbins)
    return convolve!(dest, A, β, δ)
end


function convolve!(dest :: Vector{T}, A :: BSplineEventStreamMatrix{T, S}, β :: Vector{T},  δ :: T) where {T <: Real, S <: AbstractVector{T}}
    basis = BSplineBasis(A.order, A.knots)
    nbins = ceil(Int, A.maxtime/δ)
    memlength = max(A.knots)
    membins = ceil(Int, memlength/δ)
    if length(dest) != nbins
        @error "Inexact binning: $(A.maxtime)/$δ != $(nbins)"
    end    
    for t in A.events
        for τ in 0:membins
            bs_act = bsplines(basis, τδ)
            update = sum(β[i]*bs_act[i] for i in eachindex(bs_act))
            bin_to_update = whichbin(t+τδ, δ)
            dest[bin_to_update] += update
        end
    end
    return dest
end





#Matrix(A :: EventStreamMatrix) = todense(A)
#Array(A :: EventStreamMatrix) = todense(A) 
