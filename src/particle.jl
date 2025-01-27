using StaticArrays
using PhysicalConstants.CODATA2018: c_0
using Unitful
export ParticleType, PEPlus, PEMinus, PGamma, PMuMinus, PMuPlus
export PNuE, PNuMu, PNuTau, PNuEBar, PNuMuBar, PNuTauBar, PHadronShower, PUnknown
export PLightSabre
export pdf_code, particle_shape
export Track, Cascade
export Particle, ParticleShape
export ptype_for_code
export is_neutrino
export shift_particle!

c_vac = ustrip(u"m/ns", c_0)


abstract type ParticleType end

"Subtype of ParticleType, used to classify particles. pdg_code() can be used to find the corresponding particle PDF code."
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

"Subtype of ParticleShape, used to classify signatures in detector. PEPlus, PEMinus, PGamma and PHadronShower have the shape Cascade; while PMuMinus, PMuPlus and PLightSabre have the shape Track. particle_shape() can be used to return shape of particle signature."
struct Track <: ParticleShape end
struct Cascade <: ParticleShape end
struct NoShape <: ParticleShape end


"""
    pdg_code(particle)

Return integer corresponding to the PDG code of a givern particle. ptype_for_code() can be used to find the struct of the corresponding particle.
    
# Arguments
- `particle::Type`: Particle to return PDG code of.
"""
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


"""
    ptype_for_code(code::Integer)

Return particle struct corresponding to a given PDG code. pdg_code() can be used to find the PDG code of a corresponding particle.

# Arguments
- `code::Integer`: PDG code to get particle struct of. 
"""
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
particle_shape(::Type) = NoShape()

is_neutrino(::Type) = false
is_neutrino(::Type{PNuE}) = true
is_neutrino(::Type{PNuEBar}) = true
is_neutrino(::Type{PNuMu}) = true
is_neutrino(::Type{PNuMuBar}) = true
is_neutrino(::Type{PNuTau}) = true
is_neutrino(::Type{PNuTauBar}) = true


"""
    Particle{T,PType<:ParticleType}

Struct containing information of a single particle

# Fields
- `position::SVector{3,T}`: Particle position given as a three dimensional static vector in detector coordinates
- `direction::SVector{3,T}`: Particle direction of motion given as a three dimensional static vector in detector coordinates
- `time::T`: Time of the event this particle corresponds to, given in unites of ns
- `energy::T`: Total energy of the particle given in unites of GeV
- `length::T`: Length of particle track given in unites of m
- `type::Type{PType}`: Particle type, chosen from any ParticleType subtype (PEPlus, PEMinus, PGamma, PMuPlus, PMuMinus, PNuE, PNuMu, PNuTau, PNuEBar, PNuMuBar, PNuTauBar, PHadronShower, PLightSabre, PUnknown)
"""
mutable struct Particle{T,PType<:ParticleType}
    position::SVector{3,T}
    direction::SVector{3,T}
    time::T
    energy::T
    length::T
    type::Type{PType}
end

function Particle(position::AbstractArray, direction::AbstractArray, time, energy, length, type)
    common_type = promote_type(eltype(position), eltype(direction), typeof(time), typeof(energy), typeof(length))
    return Particle(SVector{3, common_type}(position), SVector{3, common_type}(direction), common_type(time), common_type(energy), common_type(length), type)
end

function Particle(pos_x, pos_y, pos_z, dir_x, dir_y, dir_z, time, energy, length, type)
    common_type = promote_type(typeof(pos_x), typeof(pos_y), typeof(pos_z), typeof(dir_x), typeof(dir_y), typeof(dir_z), typeof(time), typeof(energy), typeof(length))
    ptype = ptype_for_code(type)
    return Particle(SVector{3, common_type}(pos_x, pos_y, pos_z), SVector{3, common_type}(dir_x, dir_y, dir_z), common_type(time), common_type(energy), common_type(length), ptype)
end



"""
    shift_particle(particle::Particle, dist_along)

Shifts the given particle along its direction by the specified distance.

# Arguments
- `particle::Particle`: The particle to be shifted.
- `dist_along`: The distance to shift the particle along its direction.

# Returns
A new `Particle` object with the shifted position.

"""
function shift_particle!(particle::Particle, dist_along)
    particle.position = particle.position + dist_along .* particle.direction
    particle.time = particle.time + dist_along / c_vac
    return particle
end

"""
    particle_shape(partricle)

Return struct corresponding to shape of argument-particle's signature
    
# Arguments
- `particle::Type`: Particle to check shape of.
"""
particle_shape(::Particle{T,PType}) where {T,PType} = particle_shape(PType)


"""
    is_neutrino(particle)

Return true if particle is of any neutrino subtype (PNuE, PNuEBar, PNuMu, PNuMuBar, PNuTau, PNuTauBar), otherwise return false.

# Arguments
- `particle::Type`: Particle to check wheather is a neutrino.
"""
is_neutrino(::Particle{T,PType}) where {T,PType} = is_neutrino(PType)


function Base.convert(::Type{Particle{T}}, x::Particle) where {T}
    pos = T.(x.position)
    dir = T.(x.direction)
    energy = T(x.energy)
    time = T(x.time)
    length = T(x.length)

    return Particle(pos, dir, time, energy, length, x.type)
end


Base.show(io::IO, ::Type{ParticleType}) = error("Invalid particle type")
Base.show(io::IO, ::Type{PMuPlus}) = print(io, "μ+")
Base.show(io::IO, ::Type{PMuMinus}) = print(io, "μ-")
Base.show(io::IO, ::Type{PEPlus}) = print(io, "e+")
Base.show(io::IO, ::Type{PEMinus}) = print(io, "e-")
Base.show(io::IO, ::Type{PGamma}) = print(io, "Photon")
Base.show(io::IO, ::Type{PNuE}) = print(io, "v_e")
Base.show(io::IO, ::Type{PNuEBar}) = print(io, "v̄_e")
Base.show(io::IO, ::Type{PNuMu}) = print(io, "v_μ")
Base.show(io::IO, ::Type{PNuMuBar}) = print(io, "v̄_μ")
Base.show(io::IO, ::Type{PNuTau}) = print(io, "v_τ")
Base.show(io::IO, ::Type{PNuTauBar}) = print(io, "v̄_τ")
Base.show(io::IO, ::Type{PUnknown}) = print(io, "Unknown Particle")

Base.show(io::IO, ::MIME"text/plain", ::Type{ParticleType}) = error("Invalid particle type")
Base.show(io::IO, ::MIME"text/plain", ::Type{PMuPlus}) = print(io, "Anti Muon")
Base.show(io::IO, ::MIME"text/plain", ::Type{PMuMinus}) = print(io, "Muon")
Base.show(io::IO, ::MIME"text/plain", ::Type{PEPlus}) = print(io, "Positron")
Base.show(io::IO, ::MIME"text/plain", ::Type{PEMinus}) = print(io, "Electron")
Base.show(io::IO, ::MIME"text/plain", ::Type{PGamma}) = print(io, "Gamma Ray")
Base.show(io::IO, ::MIME"text/plain", ::Type{PNuE}) = print(io, "Electron Neutrino")
Base.show(io::IO, ::MIME"text/plain", ::Type{PNuEBar}) = print(io, "Electron Anti-Neutrino")
Base.show(io::IO, ::MIME"text/plain", ::Type{PNuMu}) = print(io, "Muon Neutrino")
Base.show(io::IO, ::MIME"text/plain", ::Type{PNuMuBar}) = print(io, "Muon Anti-Neutrino")
Base.show(io::IO, ::MIME"text/plain", ::Type{PNuTau}) = print(io, "Tauon Neutrino")
Base.show(io::IO, ::MIME"text/plain", ::Type{PNuTauBar}) = print(io, "Tauon Anti-Neutrino")
Base.show(io::IO, ::MIME"text/plain", ::Type{PUnknown}) = print(io, "Unknown Particle")



function Base.show(io::IO, part::Particle)

    print(io, part.type, "(", "Pos:", part.position, ", Dir:", part.direction, ", T:", part.time, " ns, E:", part.energy, " GeV, L:", part.length, " m)")

end

function Base.show(io::IO, ::MIME"text/plain", part::Particle)

    show(io, MIME("text/plain"), part.type)
    print(io, " with:\n")
    print(io, "    Position: ", part.position, "\n")
    print(io, "    Direction: ", part.direction, "\n")
    print(io, "    Time: ", part.time, " ns\n")
    print(io, "    Energy: ", part.energy, " GeV\n")
    print(io, "    Length: ", part.length, " m\n")
end