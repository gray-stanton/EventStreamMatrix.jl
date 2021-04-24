
function memorylengths_away(t, reltime, memory)
    t == reltime - memory && return 0 # check exact boundary
    return ceil(Int, (nextfloat(t)-reltime)/memory)
end

function binstart(i, δ)
    i == 1 && return zero(typeof(δ))
    return (i-1)*δ
end

function binend(i, δ)
    i == 1 && return δ
    return (i)*δ
end

# Given a time and binwidth, find the time of the next bin midpoint (may be in same bin or next bin).
# if t is exactly a binmid, return t.
function nextbinmid(t, δ)
    # Need to get fancy to avoid rounding errors
    return Base._round_invstep(
        Base._round_invstep(prevfloat(t), 1/δ, RoundNearest) + δ/2,
        2/δ,
        RoundNearest
        ) 
end

function prevbinmid(t, δ)
    return Base._round_invstep(
        Base._round_invstep(nextfloat(t), 1/δ, RoundNearest) - δ/2,
        2/δ,
        RoundNearest
    )
end



function binmid(i, δ)
    i == 1 && return δ/2
    return Base._round_invstep((i-1/2)* δ, 2/δ, RoundNearest)
end

function whichbin(t :: Real, δ :: Real) 
    t == 0 && return 1
    return ceil(Int, t/δ)# Next strictly greater integer, right closed left open
end

