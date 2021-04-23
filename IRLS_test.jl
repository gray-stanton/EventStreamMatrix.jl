using EventStreamMatrix
using IterativeSolvers
using LinearAlgebra

logistic(x) = 1/(1+exp(-x))

function logistic_IRLS(E :: FirstOrderBSplineEventStreamMatrix, y :: Vector{Float64}, niter = 100)
    beta = zeros(E.ncols)
    weights = ones(E.nbins)/4
    z = y
    G = WeightedGramMatrix(E, weights)
    for i in 1:niter
        b = XtWy(E, weights, z)
        cg!(beta, G, b)
        eta = XWb(E, ones(E.ncols), beta)
        p = logistic.(eta)
        z = eta + (y - p)./weights
        weights = p .* (1 .- p)
    end
    return beta
end

function logistic_IRLS2(X :: Matrix, y :: Vector, niter = 100)
    beta = zeros(size(X)[2])
    weights = ones(size(X)[1])/4
    z = y
    G = X' * diagm(weights) * X
    for i in 1:niter
        b = X' * diagm(weights) * z
        cg!(beta, G, b)
        eta = X * beta
        p = logistic.(eta)
        z = eta + (y - p)./weights
        weights = p .* (1 .- p)
    end
    return beta
end

