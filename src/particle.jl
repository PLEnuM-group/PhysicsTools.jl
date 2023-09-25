using StaticArrays

export ParticleType, PEPlus, PEMinus, PGamma, PMuMinus, PMuPlus
export PNuE, PNuMu, PNuTau, PNuEBar, PNuMuBar, PNuTauBar, PHadronShower
export PLightSabre
export pdf_code, particle_shape
export Track, Cascade
export Particle, ParticleShape
export ptype_for_code
export is_neutrino


abstract type ParticleType end

struct PEPlus <: ParticleType end
struct PEMinus <: ParticleType end
struct PGamma <: ParticleType end
struct PMuPlus <: ParticleType end
struct PMuMinus <: ParticleType end
struct PNuE <: ParticleType end
struct PNuMu <: ParticleType end
struct PNuTau <: ParticleType end
struct PNuEBar <: ParticleType end
struct PNuMuBar <: ParticleType end
struct PNuTauBar <: ParticleType end
struct PHadronShower <: ParticleType end
struct PLightSabre <: ParticleType end
struct PUnknown <: ParticleType end

const ALL_PARTICLES = [PEPlus, PEMinus, PGamma, PMuPlus, PMuMinus, PNuE, PNuMu, PNuTau, PNuEBar, PNuMuBar, PNuTauBar, PHadronShower, PLightSabre]

abstract type ParticleShape end
struct Track <: ParticleShape end
struct Cascade <: ParticleShape end

pdg_code(::Type{PEPlus}) = -11
pdg_code(::Type{PEMinus}) = 11
pdg_code(::Type{PGamma}) = 22
pdg_code(::Type{PMuMinus}) = 13
pdg_code(::Type{PMuPlus}) = -13
pdg_code(::Type{PNuE}) = 12
pdg_code(::Type{PNuEBar}) = -12
pdg_code(::Type{PNuMu}) = 14
pdg_code(::Type{PNuMuBar}) = -14
pdg_code(::Type{PNuTau}) = 16
pdg_code(::Type{PNuTauBar}) = -16
pdg_code(::Type{PUnknown}) = 0
pdg_code(::Type) = 0

function ptype_for_code(code::Integer)
    for p in ALL_PARTICLES
        if pdg_code(p) == code
            return p
        end
    end
    return PUnknown
end

particle_shape(::Type{<:PEPlus}) = Cascade()
particle_shape(::Type{<:PEMinus}) = Cascade()
particle_shape(::Type{<:PGamma}) = Cascade()
particle_shape(::Type{<:PMuMinus}) = Track()
particle_shape(::Type{<:PMuPlus}) = Track()
particle_shape(::Type{<:PHadronShower}) = Cascade()
particle_shape(::Type{<:PLightSabre}) = Track()

is_neutrino(::Type) = false
is_neutrino(::Type{PNuE}) = true
is_neutrino(::Type{PNuEBar}) = true
is_neutrino(::Type{PNuMu}) = true
is_neutrino(::Type{PNuMuBar}) = true
is_neutrino(::Type{PNuTau}) = true
is_neutrino(::Type{PNuTauBar}) = true



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
