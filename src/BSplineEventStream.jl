struct FirstOrderBSplineEventStreamMatrix{T <: Real, L} <: AbstractEventStreamMatrix{T, L}
    events :: Vector{Tuple{T, Int}}
    labels :: Vector{L}
    δ :: T
    maxtime :: T
    nbins :: Int
    ncols :: Int
    basis :: BSplineBasis{Vector{T}}
    memory :: T
    intercept :: Bool
    function FirstOrderBSplineEventStreamMatrix(
        events,
        labels,
        δ,
        maxtime,
        splineorder,
        breakpoints,
        intercept=false
    )
        newbasis = BSplineBasis(splineorder, breakpoints)
        return new{typeof(δ), eltype(labels)}(
           events,
           labels,
           δ, 
           maxtime, 
           ceil(Int, maxtime/δ),
           length(labels) * length(newbasis) + convert(Int, intercept),
           newbasis,
           maximum(knots(newbasis)),
           intercept
    )
    end
end

#default intercept to false
FirstOrderBSplineEventStreamMatrix(events, labels, δ, maxtime, splineorder, breakpoints) = FirstOrderBSplineEventStreamMatrix(events, labels, δ, maxtime, splineorder, breakpoints, false)

function FirstOrderBSplineEventStreamMatrix(events::Vector{Tuple{Float64, String}}, labels, δ, maxtime, splineorder, breakpoints, intercept)
    relab_events = [(e, findfirst(lab -> lab==l, labels)) for (e, l) in events]
    return FirstOrderBSplineEventStreamMatrix(relab_events, labels, δ, maxtime, splineorder, breakpoints, intercept)
end

#function FirstOrderBSplineEventStreamMatrix(events::Vector{Tuple{Float64, Int}}, labels, δ, maxtime, splineorder, breakpoints, intercept)
#    return FirstOrderBSplineEventStreamMatrix{Float64, eltype(labels)}(events, labels, δ, maxtime, splineorder, breakpoints, intercept)
#end

size(E :: FirstOrderBSplineEventStreamMatrix) = (E.nbins, E.ncols)



struct FirstOrderDiscBSplineEventStreamMatrix <: AbstractEventStreamMatrix{Float64, String}
    events :: Vector{Tuple{Int, Int}}
    labels :: Vector{String}
    δ :: Float64
    maxtime :: Float64
    nbins :: Int
    ncols :: Int
    basis :: BSplineBasis{Vector{Float64}}
    memory :: Float64
    membins :: Int
    expansion :: Matrix{Float64}
    intercept :: Bool
    function FirstOrderDiscBSplineEventStreamMatrix(
        events,
        labels,
        δ,
        maxtime,
        splineorder,
        breakpoints
    )
        newbasis = BSplineBasis(splineorder, breakpoints)
        points = binstart(1, δ):δ:prevbinmid(support(newbasis)[2] + δ/2, δ)
        expmat = basismatrix(newbasis, points)
        return new(
            events,
            labels,
            δ,
            maxtime,
            ceil(Int, maxtime/δ),
            length(labels) * length(newbasis) + convert(Int, intercept),
            newbasis,
            support(newbasis)[2],
            ceil(Int, support(newbasis)[2]/δ)+1,
            expmat,
            intercept
        )
    end
end

function FirstOrderDiscBSplineEventStreamMatrix(events::Vector{Tuple{Float64, String}}, labels, δ, maxtime, splineorder, breakpoints, intercept)
    binned_events = [(whichbin(e, δ), findfirst(lab -> lab==l, labels)) for (e, l) in events]
    return FirstOrderDiscBSplineEventStreamMatrix(binned_events, labels, δ, maxtime, splineorder, breakpoints, intercept)
end

size(E :: FirstOrderDiscBSplineEventStreamMatrix) = (E.nbins, E.ncols)
setindex(E :: FirstOrderDiscBSplineEventStreamMatrix, I...) = @error "Not implemented"
function getindex(E :: FirstOrderDiscBSplineEventStreamMatrix, I...)
    bin_ind = I[1]
    if E.intercept && I[2] == 1
        return 1.0 # in intercept column
    end
    col_ind = E.intercept ? I[2] -1 : I[2]
    rellabel = E.labels[ceil(Int, col_ind/length(E.basis))]
    relspline_num = (col_ind -1) % length(E.basis) + 1
    in_range = searchsorted(E.events, bin_ind, by= (t) -> ceil(Int, (t[1]- bin_ind)/E.membins))
    out = 0.0
    for (x, l) in E.events[in_range]
        out += E.expansion[bin_ind - x + 1, relspline_num]
    end
    return out
end


struct WeightedGramMatrix{T, L} <: AbstractMatrix{T}
    X :: FirstOrderBSplineEventStreamMatrix{T, L}
    W :: Vector{T}
end

size(G :: WeightedGramMatrix) = (G.X.ncols, G.X.ncols)
Matrix(G :: WeightedGramMatrix) = XtWX(G.X, G.W)
Array(G :: WeightedGramMatrix) = Matrix(G)
setindex!(G :: WeightedGramMatrix, I...) = @error "Not implemented"
getindex(G :: WeightedGramMatrix, I...) = XtWX(G.X, G.W)[I[1], I[2]]

function mul!(ws, G::WeightedGramMatrix, b :: Vector)
    XtWXb!(ws, G.X, G.W, b)
    #ws[:] = Matrix(G.X)' * diagm(G.W) * Matrix(G.X) * b
end

function getindex(E :: FirstOrderBSplineEventStreamMatrix, I...)
    # get range of events which can influence it
    bin_ind = I[1]
    if E.intercept && I[2] == 1
        return 1.0 # in intercept column
    end
    col_ind = E.intercept ? I[2] -1 : I[2]
    rellabelnum = ceil(Int, col_ind/length(E.basis))
    reltime = binmid(bin_ind, fineness(E))
    relspline_num = (col_ind -1) % length(E.basis) + 1
    relspline = E.basis[relspline_num]
    # only spikes in range [reltime - memory, reltime] should be included
    binrange = searchsorted(events(E), prevfloat(reltime), by=(tup) -> memorylengths_away(tup[1], reltime, E.memory))
    tdiffs = [reltime - e[1] for e in events(E)[binrange] if e[2] == rellabelnum]
    out = sum(relspline.(tdiffs))
    return out
end

 Array(E :: FirstOrderDiscBSplineEventStreamMatrix) = Matrix(E)
function Matrix(E :: FirstOrderDiscBSplineEventStreamMatrix)
    out= zeros(Float64, E.nbins, E.ncols - convert(Int, E.intercept))
    nsplines = length(E.basis)
    starts = 1:nsplines:(E.ncols)
    for (x, l) in E.events
        update_bins = 1:min(E.membins, E.nbins - x )
        out[x .+ update_bins, starts[l]:(starts[l]+nsplines-1)] += E.expansion[update_bins, :]
    end
    if E.intercept
        out = hcat(ones(Float64, size(out)[1]), out)
    end
    return out
end

function setindex!(E :: FirstOrderBSplineEventStreamMatrix, I...)
    @error "Not implemented"
end

Array(E :: FirstOrderBSplineEventStreamMatrix) = Matrix(E)
function Matrix(E :: FirstOrderBSplineEventStreamMatrix) 
    #TODO make dense type adjustable
    out = zeros(Float64, E.nbins, E.ncols - convert(Int, E.intercept))
    #TODO avoid reallocating
    nsplines = length(E.basis)
    #dest = zeros(Float64, order(E.basis))
    #label_order = Dict(l => i for (i, l) in enumerate(E.labels))
    starts = 1:nsplines:(E.ncols)
    for (t, i) in events(E)
        firstpoint = nextbinmid(t, E.δ)
        lastpoint = min(prevbinmid(t + E.memory, E.δ), prevbinmid(E.maxtime, E.δ))
        points = (firstpoint:E.δ:lastpoint) .- t
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
    if E.intercept
        return hcat(ones(Float64, size(out)[1]), out)
    else
        return out
    end
end


function XWb!(dest :: Vector{T}, E :: FirstOrderBSplineEventStreamMatrix{T}, W :: Vector{T}, b :: Vector{T}) where T
    intercept = E.intercept
    #ncols = E.intercept ? E.ncols - 1 : E.ncols
    if length(dest) != E.nbins && !intercept
        throw(DimensionMismatch())
    end
    nsplines = length(E.basis)
    spl_start = E.intercept ? 2 : 1
    fo_splines = [Spline(E.basis, W[k:k+nsplines-1].*b[k:(k+nsplines-1)]) for k in spl_start:nsplines:(E.ncols)]
    bs_ws = zeros(eltype(b), order(E.basis))
    fill_val = intercept ? W[1] * b[1] : 0.0
    fill!(dest, fill_val)
    for (t, i) in events(E)
        firstpoint = nextbinmid(t, E.δ)
        lastpoint = min(prevbinmid(t + E.memory, E.δ), prevbinmid(E.maxtime, E.δ))
        points = (firstpoint:E.δ:lastpoint)
        # Broadcast over points
        for p in points
            dest[whichbin(p, E.δ)] += fo_splines[i](p - t, workspace=bs_ws)
        end
        #update_vec = fo_splines[i].(points; workspace=bs_ws)
        #update_bins = whichbin(firstpoint, E.δ):1:whichbin(lastpoint, E.δ)
        #dest[update_bins] += update_vec
    end
    return dest
end






function XWb(E :: FirstOrderBSplineEventStreamMatrix{T}, W:: Vector{T}, b :: Vector{T}) where T
    out = zeros(eltype(b), E.nbins)
    XWb!(out, E, W, b)
    return out
end

function XtWy!(dest, E :: FirstOrderBSplineEventStreamMatrix{T}, W, y :: Vector{T}) where T
    intercept = E.intercept
    if length(dest) != E.ncols && !intercept
        throw(DimensionMismatch())
    end
    nsplines = length(E.basis)
    starts = 1:nsplines:(E.ncols- convert(Int, E.intercept))
    # allocate up to maximum possible number of points
    basmat_ws = zeros(eltype(y), length(whichbin(0.0, E.δ):E.δ:whichbin(E.memory, E.δ)), nsplines)
    dest[:] .= zero(eltype(dest))
    if intercept
        dest[1] = sum(W .* y)
        dest_rest = view(dest, 2:length(dest))
    end
    for (t, i) in events(E)
        firstpoint = nextbinmid(t, E.δ)
        lastpoint = min(prevbinmid(t + E.memory, E.δ), prevbinmid(E.maxtime, E.δ))
        points = (firstpoint:E.δ:lastpoint) .- t
        basmat_view = view(basmat_ws, 1:length(points), :)
        basismatrix!(basmat_view, E.basis, points)
        update_bins = whichbin(firstpoint, E.δ):1:whichbin(lastpoint, E.δ)
        relW = W[update_bins]
        rely = y[update_bins]
        update_vec = transpose(basmat_view) * (relW .* rely)
        dest_rest[starts[i]:(starts[i] +nsplines -1)] += update_vec
    end
    return dest
end

function XtWy(E :: FirstOrderBSplineEventStreamMatrix{T}, W, y :: Vector{T}) where T
    dest = zeros(eltype(y), E.ncols)
    return XtWy!(dest, E, W, y)
end

function XtWX!(dest, E :: FirstOrderBSplineEventStreamMatrix, W)
    #stopgap no-op. Implement my bounded memory dequeue idea to avoid quadratic scaling in nspikes
    fill!(dest, zero(eltype(dest)))
    if size(dest) != (E.ncols, E.ncols)
        throw(DimensionMismatch())
    elseif length(W) != E.nbins
        throw(DimensionMismatch())
    end
    #break X into blocks of size nsplines
    nsplines = length(E.basis)
    δ = E.δ
    starts = 1:nsplines:(E.ncols- convert(Int, E.intercept))
    if E.intercept 
        dest[1, 1] = sum(W)
        dest_rest = view(dest, 2:size(dest)[1], 2:size(dest)[2])
        # should have first row and first column be sum of w * basis values in those columns...
    else
        dest_rest = dest
    end
    exp_mats = Queue{Tuple{Float64, Int, Matrix{Float64}}}()
    curevent_idx = 0
    maxmemevent_idx = 0
    while curevent_idx < length(E.events)
        #Queue setup
        curevent_idx += 1
        # No more events to queue
        newt, _ = E.events[curevent_idx]
        for (t, j) in E.events[(maxmemevent_idx+1):end]
            if t - newt < E.memory
                firstpoint = nextbinmid(t, E.δ)
                lastpoint = min(prevbinmid(t + E.memory, E.δ), prevbinmid(E.maxtime, E.δ))
                enqueue!(exp_mats, (t, j, basismatrix(E.basis, (firstpoint:δ:lastpoint) .- t)))
                maxmemevent_idx += 1
            else
                break
            end
        end

        # Handle self Mult
        t1, j1, fullb1 = dequeue!(exp_mats)
        memstart_point = nextbinmid(t1, E.δ)
        memend_point = min(prevbinmid(t1 + E.memory, E.δ), prevbinmid(E.maxtime, E.δ))
        update_bins = whichbin(memstart_point, E.δ):1:whichbin(memend_point, E.δ)
        b1 = view(fullb1, 1:length(update_bins), :)
        update_cols1 = starts[j1]:(starts[j1]+nsplines - 1)
        relW = Diagonal(W[update_bins])
        dest_rest[update_cols1, update_cols1] += b1' * relW * b1
        # Add intercept 
        if E.intercept
            intsum = sum(b1 .* W[update_bins]; dims=1)
            dest[[1], update_cols1 .+ 1] += reshape(intsum, (1, length(update_cols1)))
            dest[update_cols1 .+ 1, [1]] += reshape(intsum, (length(update_cols1), 1))
        end
        # Handle interaction with all other events in memory
        for (t2, j2, b2) in exp_mats
            interactstart_point = nextbinmid(t2, E.δ)
            update_bins = whichbin(interactstart_point, E.δ):1:whichbin(memend_point, E.δ)
            update_cols2 = starts[j2]:(starts[j2]+nsplines-1)
            relW = Diagonal(W[update_bins])
            relb1 = view(b1, (size(b1)[1]-length(update_bins)+1):size(b1)[1],:)
            relb2 = view(b2, 1:length(update_bins), :)
            update_mat1 = relb1' * relW * relb2
            update_mat2 = transpose(update_mat1)
            if j1 == j2
                # Diagonal block of output
                dest_rest[update_cols1, update_cols1] += update_mat1
                dest_rest[update_cols1, update_cols1] += update_mat2
            else
                # Off-diagonal block of output
                dest_rest[update_cols1, update_cols2] += update_mat1
                dest_rest[update_cols2, update_cols1] += update_mat2
            end
        end
        # Expand new matrices, add to queue
    end
    return dest
end



function XtWX!_old(dest, E, W)
    if size(dest) != (E.ncols, E.ncols)
        throw(DimensionMismatch())
    elseif length(W) != E.nbins
        throw(DimensionMismatch())
    end
    label_order = Dict(l => i for (i, l) in enumerate(E.labels))
    nsplines = length(E.basis)
    starts = 1:nsplines:(E.ncols)
    basmat_ws1 = zeros(eltype(W), length(whichbin(0.0, E.δ):E.δ:whichbin(E.memory, E.δ)), nsplines)
    basmat_ws2 = zeros(eltype(W), length(whichbin(0.0, E.δ):E.δ:whichbin(E.memory, E.δ)), nsplines)
    dest[:] .= zero(eltype(dest)) # Zero out matrix
    for (j, (t1, l1)) in enumerate(events(E))
        self_mult = true # first case is always self, handle differently.
        for (t2, l2) in events(E)[j:end]
            if t2 > t1 + E.memory
                break
            else
                i1 = label_order[l1]
                i2 = label_order[l2]
                firstpoint = max(nextbinmid(t1, E.δ), nextbinmid(t2, E.δ))
                lastpoint = min(prevbinmid(t1 + E.memory, E.δ),  prevbinmid(E.maxtime, E.δ))
                points_for_t1 = (firstpoint:E.δ:lastpoint) .- t1
                points_for_t2 = (firstpoint:E.δ:lastpoint) .- t2
                bmat_for_t1 = view(basmat_ws1, 1:length(points_for_t1), :)
                bmat_for_t2 = view(basmat_ws2, 1:length(points_for_t2), :)
                basismatrix!(bmat_for_t1, E.basis, points_for_t1)
                basismatrix!(bmat_for_t2, E.basis, points_for_t2)
                update_bins = whichbin(firstpoint, E.δ):1:whichbin(lastpoint, E.δ)
                relW = W[update_bins]
                # This is the "above" diagonal block, hence want the smaller of the two to transpose
                update_mat1 = bmat_for_t1' * Diagonal(relW) * bmat_for_t2 # TODO: Make this into a view as well.
                update_mat2 = transpose(update_mat1)
                # Symmetric, update both blocks
                # If same event, only update once. Otherwise there are two contributions.
                if i1 == i2
                    dest[starts[i1]:(starts[i1]+nsplines-1), starts[i2]:(starts[i2] + nsplines-1)] += update_mat1
                    if !self_mult
                        dest[starts[i1]:(starts[i1]+nsplines-1), starts[i1]:(starts[i1] + nsplines-1)] += update_mat2
                    end
                else
                    #k1 = min(i1, i2)
                    #k2 = max(i1, i2)
                    dest[starts[i1]:(starts[i1]+nsplines-1), starts[i2]:(starts[i2] + nsplines-1)] += update_mat1
                    dest[starts[i2]:(starts[i2]+nsplines-1), starts[i1]:(starts[i1] + nsplines-1)] += update_mat2
                end
                self_mult = false 
            end
        end
    end
    return dest
end


#function XtWX(E:: FirstOrderBSplineEventStreamMatrix, W=ones(E.nbins))
#    out = zeros(Float64, E.ncols, E.ncols)
#    XtWX!(out, E, W)
#    return out
#end

function XtWXb!(dest, E :: FirstOrderBSplineEventStreamMatrix, W, b ::Vector{T}) where T
    # stopgap implementation as XtW(XB)
    if length(dest) != size(E)[2]
        throw(DimensionMismatch())
    elseif length(b) != size(E)[2]
        throw(DimensionMismatch())
    elseif length(W) != size(E)[1]
        throw(DimensionMismatch())
    else
        xb = XWb(E, ones(size(E)[2]), b)
        XtWy!(dest, E, W, xb)
        return dest
    end
end

function XtWXb!_old(dest , E::FirstOrderBSplineEventStreamMatrix, W, b::Vector{T}) where T
    if length(dest) != E.ncols
        throw(DimensionMismatch())
    elseif length(W) != E.nbins
        throw(DimensionMismatch())
    end
    label_order = Dict(l => i for (i, l) in enumerate(E.labels))
    nsplines = length(E.basis)
    starts = 1:nsplines:(E.ncols)
    dest[:] .= zero(eltype(dest)) 
    for (j, (t1, l1)) in enumerate(events(E))
        self_mult = true
        for (t2, l2) in events(E)[j:end]
            if t2 > t1 + E.memory
                break
            else
                i1 = label_order[l1]
                i2 = label_order[l2]
                firstpoint = max(nextbinmid(t1, E.δ), nextbinmid(t2, E.δ))
                lastpoint = min(prevbinmid(t1 + E.memory, E.δ),  prevbinmid(E.maxtime, E.δ))
                points_for_t1 = (firstpoint:E.δ:lastpoint) .- t1
                points_for_t2 = (firstpoint:E.δ:lastpoint) .- t2
                bmat_for_t1 = basismatrix(E.basis, points_for_t1)
                bmat_for_t2 = basismatrix(E.basis, points_for_t2)
                update_bins = whichbin(firstpoint, E.δ):1:whichbin(lastpoint, E.δ)
                relW = W[update_bins]
                relb1 = b[starts[i1]:(starts[i1] + nsplines -1)]
                relb2 = b[starts[i2]:(starts[i2] + nsplines -1)]
                B2beta = bmat_for_t2 * relb2
                B1beta = bmat_for_t1 * relb1
                #update_mat1 = bmat_for_t1' * Diagonal(relW) * bmat_for_t2 
                #update_mat2 = transpose(update_mat1)
                #update_vec1 = update_mat * b[starts[i1]:(starts[i1] + nsplines -1)]
                #update_vec2 = transpose(update_mat) * b[starts[i2]:(starts[i2] + nsplines -1)]
                #beta_up1 = update_mat1 * relb2
                beta_up1 = bmat_for_t1'*Diagonal(relW)*B2beta
                beta_up2 = bmat_for_t2'*Diagonal(relW)*B1beta
                #beta_up2 = update_mat2 * relb1
                if i1 == i2
                    dest[starts[i1]:(starts[i1] + nsplines -1)] += beta_up1
                    if !self_mult
                        dest[starts[i1]:(starts[i1] +nsplines -1)] += beta_up2
                    end
                else
                    # Avoid double-counting 
                    dest[starts[i1]:(starts[i1] + nsplines -1)] += beta_up1
                    dest[starts[i2]:(starts[i2] + nsplines -1)] += beta_up2
                end
                self_mult = false
            end
        end
    end
    return dest
end




function XtWXb(E :: FirstOrderBSplineEventStreamMatrix{T}, W, b :: Vector{T}) where T
    out = zeros(Float64, E.ncols)
    XtWXb!(out, E, W, b)
    return out
end





#= function convolve(A :: BSplineEventStreamMatrix{T, S}, β :: Vector{T},  δ :: T) where {T <: Real, S <: AbstractVector{T}}
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
 =#




#Matrix(A :: EventStreamMatrix) = todense(A)
#Array(A :: EventStreamMatrix) = todense(A) 
