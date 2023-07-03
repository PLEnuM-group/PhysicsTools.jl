using StaticArrays

export ParticleType, PEPlus, PEMinus, PGamma, PMuMinus, PMuPlus
export pdf_code, particle_shape
export Track, Cascade
export Particle, ParticleShape


abstract type ParticleType end

struct PEPlus <: ParticleType end
struct PEMinus <: ParticleType end
struct PGamma <: ParticleType end
struct PMuPlus <: ParticleType end
struct PMuMinus <: ParticleType end

abstract type ParticleShape end
struct Track <: ParticleShape end
struct Cascade <: ParticleShape end

pdg_code(::Type{PEPlus}) = -11
pdg_code(::Type{PEMinus}) = 11
pdg_code(::Type{PGamma}) = 22
pdg_code(::Type{PMuMinus}) = 13
pdg_code(::Type{PMuPlus}) = -13

particle_shape(::Type{<:PEPlus}) = Cascade()
particle_shape(::Type{<:PEMinus}) = Cascade()
particle_shape(::Type{<:PGamma}) = Cascade()
particle_shape(::Type{<:PMuMinus}) = Track()
particle_shape(::Type{<:PMuPlus}) = Track()

mutable struct Particle{T,PType<:ParticleType}
    position::SVector{3,T}
    direction::SVector{3,T}
    time::T
    energy::T
    length::T
    type::Type{PType}
end

particle_shape(::Particle{T,PType}) where {T, PType} = particle_shape(PType)


function Base.convert(::Type{Particle{T}}, x::Particle) where {T}
    pos = T.(x.position)
    dir = T.(x.direction)
    energy = T(x.energy)
    time = T(x.time)
    length = T(x.length)

    return Particle(pos, dir, time, energy, length, x.type)
end
