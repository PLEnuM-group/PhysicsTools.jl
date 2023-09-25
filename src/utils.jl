module Utils
using StaticArrays
using FastGaussQuadrature
using LinearAlgebra
using DataStructures
using Distributions
using Roots
using Polynomials
using Random
using PhysicalConstants.CODATA2018

export fast_linear_interp, transform_integral_range
export integrate_gauss_quad
export sph_to_cart, apply_rot, cart_to_sph, rot_to_ez_fast, rot_from_ez_fast, calc_rot_matrix
export CategoricalSetDistribution
export sample_cherenkov_track_direction
export rand_gamma
export fwhm, calc_gamma_shape_mean_fwhm
export repeat_for, repeat_for!, split_by
export ssc
export gumbel_width_from_fwhm
export cart_to_cyl, cart_to_cyl
export frank_tamm, frank_tamm_norm

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
    fast_linear_interp(x_eval::Number, xs::AbstractVector, ys::AbstractVector)

Linearly interpolate xs -> ys and evaluate x_eval on interpolation. Assume xs are sorted in ascending order.
"""
function fast_linear_interp(x_eval::Number, xs::AbstractVector, ys::AbstractVector)

    if length(xs) != length(ys)
        error("Input vectors must have the same length")
    end

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

"""
    fast_linear_interp(x_eval::Number, knots::AbstractVector, lower::Number, upper::Number)

Linearly interpolate knots and evaluate x_eval on interpolation.
Assume knots are equidistant in (lower, upper).
"""
function fast_linear_interp(x_eval::Number, knots::AbstractVector, lower::Number, upper::Number)

    x_eval = clamp(x_eval, lower, upper)
    range = upper - lower
    n_knots = size(knots, 1)
    step_size = range / (n_knots - 1)

    along_range = (x_eval - lower) / step_size
    along_range_floor = floor(along_range)
    lower_knot = Int64(along_range_floor) + 1

    if lower_knot == n_knots
        return @inbounds oftype(along_range, knots[end])
    end

    along_step = along_range - along_range_floor
    @inbounds y_low = knots[lower_knot]
    @inbounds slope = (knots[lower_knot+1] - y_low)

    interpolated = y_low + slope * along_step

    return interpolated
end

"""
    transform_integral_range(x::Real, f::T, xrange::Tuple{<:Real,<:Real}) where {T<:Function}

Apply change of interval formula for evaluating quadrature rule on `f` in interval `xrange`.
See: https://en.wikipedia.org/wiki/Gaussian_quadrature#Change_of_interval
"""
function transform_integral_range(x::Real, f::T, xrange::Tuple{<:Real,<:Real}) where {T<:Function}
    ba_half = (xrange[2] - xrange[1]) / 2
    x = oftype(ba_half, x)

    u_traf = ba_half * x + (xrange[1] + xrange[2]) / 2
    return oftype(x, f(u_traf) * ba_half)
end

"""
    integrate_gauss_quad(f::T, a::Real, b::Real) where {T<:Function}
Integrate `f` in interval [`a`, `b`] using 10-point Gauss-Legendre quadrature.
"""
function integrate_gauss_quad(f::T, a::Real, b::Real) where {T<:Function}
    U = promote_type(typeof(a), typeof(b))
    return U(integrate_gauss_quad(f, a, b, GL10[1], GL10[2]))
end

"""
    integrate_gauss_quad(f::T, a::Real, b::Real, order::Integer) where {T<:Function}
Integrate `f` in interval [`a`, `b`] using Gauss-Legendre quadrature of order `order`.
"""
function integrate_gauss_quad(f::T, a::Real, b::Real, order::Integer) where {T<:Function}
    nodes, weights = gausslegendre(order)
    integrate_gauss_quad(f, a, b, nodes, weights)
end

"""
    integrate_gauss_quad(f::T, a::Real, b::Real, nodes::AbstractVector{U}, weights::AbstractVector{U}) where {T<:Function,U<:Real}
Integrate `f` in interval [`a`, `b`] using Gaussian quadrature with `nodes` and `weights`.
"""
function integrate_gauss_quad(f::T, a::Real, b::Real, nodes::AbstractVector{U}, weights::AbstractVector{U}) where {T<:Function,U<:Real}    
    dot(weights, map(x -> transform_integral_range(x, f, (a, b)), nodes))
end

"""
    sph_to_cart(theta::Real, phi::Real)
Convert spherical to cartesian coordinates.

Uses ISO convention (inclination, azimuth).
"""
function sph_to_cart(theta::Real, phi::Real)

    sin_theta, cos_theta = sincos(theta)
    sin_phi, cos_phi = sincos(phi)

    T = promote_type(typeof(theta), typeof(phi))
    x::T = cos_phi * sin_theta
    y::T = sin_phi * sin_theta
    z::T = cos_theta

    return SA[x, y, z]
end

sph_to_cart(x::AbstractVector) = sph_to_cart(x[1], x[2])

"""
    cart_to_sph(x::Real, y::Real, z::Real)

Convert cartesian to spherical coordinates. Assumes x, y, z represent a unit vector.
Uses ISO convetion (inclination, azimuth).
"""
function cart_to_sph(x::Real, y::Real, z::Real)

    T = promote_type(typeof(x), typeof(y), typeof(z))
    z = clamp(z, -1, 1)
    if z == 1
        return SA{T}[0, 0]
    elseif z == -1
        return SA{T}[π, 0]
    end

    theta = acos(z)    
    phi = atan(y, x)

    phi = phi < 0 ? phi + 2 * π : phi
    phi = phi > 2 * π ? phi - 2 * π : phi

    return SA{T}[theta, phi]
end

cart_to_sph(x::AbstractVector) = cart_to_sph(x[1], x[2], x[3])

"""
    cart_to_cyl(x::Real, y::Real, z::Real)
Convert cartesian to cylinder coordinates (rho, phi, z)
"""
function cart_to_cyl(x::Real, y::Real, z::Real)

    T = promote_type(typeof(x), typeof(y), typeof(z))

    rho = sqrt(x^2 + y^2)
    phi = acos(x / rho)
    phi = y >= 0 ? phi : 2*π - phi
    return SA{T}[rho, phi, z]
end

cart_to_cyl(x::AbstractArray) = cart_to_cyl(x[1], x[2], x[3])

"""
    cyl_to_cart(rho::Real, phi::Real, z::Real)
Convert cylinder (rho, phi, z) to cartesian coordinated
"""
function cyl_to_cart(rho::Real, phi::Real, z::Real)
    T = promote_type(typeof(rho), typeof(phi), typeof(z))
    return SA{T}[rho * cos(phi), rho*sin(phi), z]
end

cyl_to_cart(x::AbstractArray) = cyl_to_cart(x[1], x[2], x[3])



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

function Base.:(==)(a::CategoricalSetDistribution, b::CategoricalSetDistribution)
    return (collect(a.set) == collect(b.set)) && (a.cat == b.cat)
end

Base.rand(rng::AbstractRNG, pdist::CategoricalSetDistribution) = pdist.set[rand(rng, pdist.cat)]
Base.rand(pdist::CategoricalSetDistribution) = rand(default_rng(), pdist)


ssc(v::AbstractVector) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

"""
    calc_rot_matrix(a, b)
Calculates rotation matrix obtained by rotating `a` to `b`.
"""
function calc_rot_matrix(a, b)
    # Rotate a to b and apply to operand
    I3 = SMatrix{3,3}(I)
    if a == b
        return I3
    end

    cross_ab = cross(a, b)
    ssc_cross_ab = ssc(cross_ab)

    R = SMatrix{3,3}(I) + ssc_cross_ab + ssc_cross_ab^2 * (1 - dot(a, b)) / norm(cross_ab)^2
    return R
end


"""
    apply_rot(a, b, operand)
Calculates rotation matrix obtained by rotating `a` to `b`. Applies resulting matrix to `operand``.
"""
function apply_rot(a, b, operand)
    R = calc_rot_matrix(a, b)
    res = R * operand
    return res ./ norm(res)
end


"""
    rot_to_ez_fast(a::AbstractVector{T}, operand::AbstractVector{T}) where {T<:Real}

Calc rotation matrix which rotates `a` to e_z. Applies resulting matrix to `operand`.
"""
@inline function rot_to_ez_fast(a::AbstractVector{T}, operand::AbstractVector{T}) where {T<:Real}

    if abs(a[3]) == T(1)
        return @SVector[operand[1], copysign(operand[2], a[3]), copysign(operand[3], a[3])]
    end

    a1sq = a[1]^2
    a2sq = a[2]^2
    fact = (1 - a[3]) / (a2sq + a1sq)
    a1a2 = a[1] * a[2]

    x = fma(-a1sq, fact, 1) * operand[1] - a1a2 * fact * operand[2] - a[1] * operand[3]
    y = -a1a2 * fact * operand[1] + fma(-a2sq, fact, 1) * operand[2] - a[2] * operand[3]
    z = a[1] * operand[1] + a[2] * operand[2] + a[3] * operand[3]
    return SA[x, y, z]
end

"""
    rot_to_ez_fast(a::SVector{3,T}, operand::SVector{3,T}) where {T<:Real}

Calc rotation matrix which rotates e_z to `a`. Applies resulting matrix to `operand`.
"""
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
    sample_cherenkov_track_direction(T::Type)

Sample direction of a Cherenkov track in a particle cascade.

Integrated over an entire particle cascade, this approximates the emission direction of Cherenkov photons.
Taken from: https://arxiv.org/abs/1301.5361, p29
"""
function sample_cherenkov_track_direction(T::Type)
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
    c::T = one(T) / sqrt(9 * d)


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

"""
    calc_gamma_shape_mean_fwhm(mean, target_fwhm)

Calculate distribution parameters `alpha` and `theta` for a Gamma distribution from
desired FWHM
"""
function calc_gamma_shape_mean_fwhm(mean, target_fwhm)
    function _optim(theta)
        alpha = mean / theta
        tt_dist = Gamma(alpha, theta)
        fwhm(tt_dist, mode(tt_dist); xlims=(0, 100)) - target_fwhm
    end

    find_zero(_optim, [0.1 * target_fwhm^2 / mean, 10 * target_fwhm^2 / mean], A42())
end


"""
    fit_gumbel_fwhm_width()

Fit a polynomial to the relationship between Gumbel width and FWHM
"""
function fit_gumbel_fwhm_width()
    # find relationship between Gumbel width and FWHM
    widths = 0.5:0.01:5
    # Fit the function width = a * fwhm + b
    poly = Polynomials.fit(map(w -> fwhm(Gumbel(0, w), w), widths), widths, 1)
    poly
end

"""
    gumbel_width_from_fwhm(theta)
Return FWHM of a Gumbel distribution with parameters `mu`=0, `theta`
"""
gumbel_width_from_fwhm = fit_gumbel_fwhm_width()


"""
    frank_tamm(wavelength::Real, ref_index::T) where {T<:Real}

Evaluate Frank-Tamm formula
"""
function frank_tamm(wavelength::Real, ref_index::T) where {T<:Real}
    return T(2 * pi * FineStructureConstant / wavelength^2 * (1 - 1 / ref_index^2))
end

"""
    frank_tamm_norm(wl_range::Tuple{T, T}, ref_index_func::Function) where {T<:Real}

Calculate number of Cherenkov photons per length in interval `wl_range`.
Returned number is in units m^-1.
"""
function frank_tamm_norm(wl_range::Tuple{T,T}, ref_index_func::Function) where {T<:Real}
    f(x) = frank_tamm(x, ref_index_func(x))
    integrate_gauss_quad(f, wl_range[1], wl_range[2]) * T(1E9)
end



end
