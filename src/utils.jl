module Utils
using StaticArrays
using FastGaussQuadrature
using LinearAlgebra
using DataStructures
using Distributions
using Roots

export fast_linear_interp, transform_integral_range
export integrate_gauss_quad
export sph_to_cart, apply_rot, cart_to_sph, rot_to_ez_fast, rot_from_ez_fast, calc_rot_matrix
export CategoricalSetDistribution
export sample_cherenkov_track_direction
export rand_gamma
export fwhm
export repeat_for, repeat_for!, split_by
export ssc

const GL10 = gausslegendre(10)

"""
    repeat_for(x::AbstractMatrix, n::AbstractVector)
Repeat each slice along the second dimension of x for n times
"""
repeat_for(x::AbstractMatrix, n::AbstractVector) = repeat_for!(x, n, similar(x, (size(x, 1), sum(n))))

function repeat_for!(x::AbstractMatrix, n::AbstractVector, out)
    ix = firstindex(x, 2)
    for i in eachindex(n)
        ni = n[i]
        out[:, ix:ix+(ni-1)] .= x[:, i]
        ix += ni
    end
    return out
end


"""
    split_by(x::AbstractVector, n::AbstractVector)
Split vector x into parts, where the split indices are given by vector n.

"""
function split_by(x::AbstractVector, n::AbstractVector)
    result = Vector{Vector{eltype(x)}}()
    start = firstindex(x)
    for l in n
        push!(result, x[start:start+l-1])
        start += l
    end

    return result
end

"""
    split_by!(x::AbstractVector, n::AbstractVector, out)
Write output into out
"""
function split_by!(x::AbstractVector, n::AbstractVector, out)
    start = firstindex(x)
    for (i, l) in enumerate(n)
        out[i] = x[start:start+l-1]
        start += l
    end
end


"""
    fast_linear_interp(x_eval::T, xs::AbstractVector{T}, ys::AbstractVector{T})

Linearly interpolate xs -> ys and evaluate x_eval on interpolation. Assume xs are sorted in ascending order.
"""
function fast_linear_interp(x_eval::T, xs::AbstractVector{T}, ys::AbstractVector{T}) where {T}

    lower = first(xs)
    upper = last(xs)
    x_eval = clamp(x_eval, lower, upper)

    if x_eval == lower
        return ys[1]
    elseif x_eval == upper
        return ys[end]
    end

    ix_upper = searchsortedfirst(xs, x_eval)
    ix_lower = ix_upper - 1

    @inbounds edge_l = xs[ix_lower]
    @inbounds edge_h = xs[ix_upper]

    step = edge_h - edge_l

    along_step = (x_eval - edge_l) / step

    @inbounds y_low = ys[ix_lower]
    @inbounds slope = (ys[ix_upper] - y_low)

    interpolated = y_low + slope * along_step

    return interpolated

end


function fast_linear_interp(x::T, knots::AbstractVector{T}, lower::T, upper::T) where {T}
    # assume equidistant

    x = clamp(x, lower, upper)
    range = upper - lower
    n_knots = size(knots, 1)
    step_size = range / (n_knots - 1)

    along_range = (x - lower) / step_size
    along_range_floor = floor(along_range)
    lower_knot = Int64(along_range_floor) + 1

    if lower_knot == n_knots
        return @inbounds knots[end]
    end

    along_step = along_range - along_range_floor
    @inbounds y_low = knots[lower_knot]
    @inbounds slope = (knots[lower_knot+1] - y_low)

    interpolated = y_low + slope * along_step

    return interpolated
end


function transform_integral_range(x::Real, f::T, xrange::Tuple{<:Real,<:Real}) where {T<:Function}
    ba_half = (xrange[2] - xrange[1]) / 2
    x = oftype(ba_half, x)

    u_traf = ba_half * x + (xrange[1] + xrange[2]) / 2
    oftype(x, f(u_traf) * ba_half)

end

function integrate_gauss_quad(f::T, a::Real, b::Real) where {T<:Function}
    U = promote_type(typeof(a), typeof(b))
    U(integrate_gauss_quad(f, a, b, GL10[1], GL10[2]))
end

function integrate_gauss_quad(f::T, a::Real, b::Real, order::Integer) where {T<:Function}
    nodes, weights = gausslegendre(order)
    integrate_gauss_quad(f, a, b, nodes, weights)
end

function integrate_gauss_quad(f::T, a::Real, b::Real, nodes::AbstractVector{U}, weights::AbstractVector{U}) where {T<:Function,U<:Real}
    dot(weights, map(x -> transform_integral_range(x, f, (a, b)), nodes))
end

function sph_to_cart(theta::Real, phi::Real)
    sin_theta, cos_theta = sincos(theta)
    sin_phi, cos_phi = sincos(phi)

    T = promote_type(typeof(theta), typeof(phi))
    x::T = cos_phi * sin_theta
    y::T = sin_phi * sin_theta
    z::T = cos_theta

    return SA[x, y, z]
end


function cart_to_sph(x::Real, y::Real, z::Real)

    T = promote_type(typeof(x), typeof(y), typeof(z))
    if z == 1
        return zero(T), zero(T)
    elseif z == -1
        return T(π), zero(T)
    end
    theta = acos(z)
    if (x == 0) && (y > 0)
        phi = T(π / 2)
    elseif (x == 0) && (y < 0)
        phi = T(-π / 2)
    else
        phi = atan(y, x)
    end

    if phi < 0
        phi = 2 * π + phi
    end

    return theta, phi
end

cart_to_sph(x::SVector{3,<:Real}) = cart_to_sph(x[1], x[2], x[3])



"""
CategoricalSetDistribution{T, U<:Real}

Represents a Categorical distribution on a set

### Examples

- `p = CategoricalSetDistribution(Set([:EMinus, :EPlus]), Categorical([0.1, 0.9]))
   rand(p)` -- returns `:EMinus` with 10% probability and `:Eplus` with 90% probability

- `p = CategoricalSetDistribution(Set([:EMinus, :EPlus]), [0.1, 0.9])` -- convenience constructor
"""
struct CategoricalSetDistribution{T}
    set::OrderedSet{T}
    cat::Categorical

    function CategoricalSetDistribution(set::OrderedSet{T}, cat::Categorical) where {T}
        if length(set) != ncategories(cat)
            error("Set and categorical have to be of same length")
        end
        new{T}(set, cat)
    end

    function CategoricalSetDistribution(set::OrderedSet{T}, probs::Vector{<:Real}) where {T}
        new{T}(set, Categorical(probs))
    end
end

Base.rand(pdist::CategoricalSetDistribution) = pdist.set[rand(pdist.cat)]


ssc(v::AbstractVector) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]


function calc_rot_matrix(a, b)
    # Rotate a to b and apply to operand
    I3 = SMatrix{3,3}(I)
    if a == b
        return I3
    end

    cross_ab = cross(a, b)
    ssc_cross_ab = ssc(cross_ab)

    R = SMatrix{3,3}(I) + ssc_cross_ab + ssc_cross_ab^2 * (1 - dot(a, b)) / norm(cross_ab)^2
end


"""
    apply_rot(a, b, operand)

Calculates rotation matrix obtained by rotating a to b. Apply to operand.
Apply the resulting rotation to operand.
"""
function apply_rot(a, b, operand)
    R = calc_rot_matrix(a, b)
    res = R * operand
    return res ./ norm(res)
end


@inline function rot_to_ez_fast(a::SVector{3,T}, operand::SVector{3,T}) where {T<:Real}

    if abs(a[3]) == T(1)
        return @SVector[operand[1], copysign(operand[2], a[3]), copysign(operand[3], a[3])]
    end

    #=
    cross_ab = [a[2] -a[1] 0]
    norm_cross_ab_sq = a[2]^2 + a[1]^2
    ssc_ab = [0 0 -a[1]; 0 0 -a[2]; a[1] a[2] 0 ]
    ssc_ab_sq = [-a[1]^2 -a[1]*a[2] 0; -a[1]*a[2] -a[2]^2 0; 0 0 -a[1]^2 -a[2]^2]

    R = [1-a1sq*fact  -a[1]*a[2]*fact -a[1] ;
        -a[1]*a[2]*fact 1-fact*a2sq   -a[2] ;
         a[1]           a[2]             a[3]]
    =#
    a1sq = a[1]^2
    a2sq = a[2]^2
    fact = (1 - a[3]) / (a2sq + a1sq)
    a1a2 = a[1] * a[2]

    x = fma(-a1sq, fact, 1) * operand[1] - a1a2 * fact * operand[2] - a[1] * operand[3]
    y = -a1a2 * fact * operand[1] + fma(-a2sq, fact, 1) * operand[2] - a[2] * operand[3]
    z = a[1] * operand[1] + a[2] * operand[2] + a[3] * operand[3]
    return SA[x, y, z]
end


@inline function rot_from_ez_fast(a::SVector{3,T}, operand::SVector{3,T}) where {T<:Real}

    if abs(a[3]) == T(1)
        return @SVector[operand[1], copysign(operand[2], a[3]), copysign(operand[3], a[3])]
    end

    a1sq = a[1]^2
    a2sq = a[2]^2
    fact = (1 - a[3]) / (a2sq + a1sq)
    a1a2 = a[1] * a[2]

    x = fma(-a1sq, fact, 1) * operand[1] - a1a2 * fact * operand[2] + a[1] * operand[3]
    y = -a1a2 * fact * operand[1] + fma(-a2sq, fact, 1) * operand[2] + a[2] * operand[3]
    z = -a[1] * operand[1] - a[2] * operand[2] + a[3] * operand[3]
    return SA[x, y, z]
end



"""
"""
function sample_cherenkov_track_direction(T::Type)
    # Mystery values from clsim
    angularDist_a = T(0.39)
    angularDist_b = T(2.61)
    angularDist_I = T(1) - exp(-angularDist_b * 2^angularDist_a)

    costheta = max(T(1) - (-log(T(1) - rand(T) * angularDist_I) / angularDist_b)^(1 / angularDist_a), T(-1))
    phi = T(2 * π) * rand(T)

    return sph_to_cart(acos(costheta), phi)

end

"""
    rand_gamma(shape, scale)

Sample gamma variates when shape > 1
"""
function rand_gamma(shape::Real, scale::Real, T::Type{U}=Float64) where {U<:Real}

    d = T(shape - 1 / 3)
    c = one(T) / sqrt(9 * d)


    while true
        x = randn(T)
        v = fma(c, x, one(T))

        if v <= 0
            continue
        end

        v = v * v * v
        u = rand(T)

        xsq = x * x

        if u < one(T) - fma(T(0.0331), xsq * xsq, one(xsq)) # 1 - 0.0331 * xsq^2
            return d * v * scale
        end

        if log(u) < T(0.5) * xsq + d * (one(T) - v + log(v))
            return d * v * scale
        end
    end
end

"""
    fwhm(d::UnivariateDistribution, xmode::Real; xlims=(-20, 20))

Calculate FWHM of a univariate distribution
"""
function fwhm(d::UnivariateDistribution, xmode::Real; xlims=(-20, 20))
    ymode = pdf(d, xmode)

    z0 = find_zero(x -> pdf(d, x) - ymode / 2, (xlims[1], xmode), A42())
    z1 = find_zero(x -> pdf(d, x) - ymode / 2, (xmode, xlims[2]), A42())
    return z1 - z0
end


end
