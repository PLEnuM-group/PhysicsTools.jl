module ProposalInterface
using PythonCall
using StaticArrays
using Logging
import Pkg

import ..Particle, ..PEMinus, ..PEPlus, ..PMuMinus, ..PMuPlus, ..PHadronShower, ..ParticleType


export proposal_secondary_to_particle, propagate_muon

const pp = PythonCall.pynew()

const PKGDIR = pkgdir(@__MODULE__)

function __init__()
    
    tmp = pyimport("proposal")
   
    if !isdir(joinpath(PKGDIR, "assets/proposal_tables"))
        mkpath(joinpath(PKGDIR, "assets/proposal_tables"))
    end
    tmp.InterpolationSettings.tables_path = joinpath(PKGDIR, "assets/proposal_tables")
    PythonCall.pycopy!(pp, tmp)
end


function proposal_secondary_to_particle(loss)
    energy = pyconvert(Float64, loss.energy) / 1E3
    pos = SA[
        pyconvert(Float64, loss.position.x/100),
        pyconvert(Float64, loss.position.y/100),
        pyconvert(Float64, loss.position.z/100)
        ]
    dir = SA[
        pyconvert(Float64, loss.direction.x),
        pyconvert(Float64, loss.direction.y),
        pyconvert(Float64, loss.direction.z)]
    time = pyconvert(Float64, loss.time) * 1E9

    pname = pyconvert(String, pp.particle.Interaction_Type(loss.type).name)

    if pname == "photonuclear"
        ptype = PHadronShower
    else
        ptype = PEMinus
    end

    return Particle(pos, dir, time, energy, 0.0, ptype)
end

function proposal_continuous_to_particle(loss)
    energy = pyconvert(Float64, loss.energy) / 1E3
    pos = SA[
        pyconvert(Float64, loss.start_position.x/100),
        pyconvert(Float64, loss.start_position.y/100),
        pyconvert(Float64, loss.start_position.z/100)
        ]
    dir = SA[
        pyconvert(Float64, loss.direction_initial.x),
        pyconvert(Float64, loss.direction_initial.y),
        pyconvert(Float64, loss.direction_initial.z)]
    time = pyconvert(Float64, loss.time_initial) * 1E9

    ptype = PEPlus

    return Particle(pos, dir, time, energy, 0.0, ptype)
end


function make_propagator(ptype::Type{<:ParticleType})
    if ptype == PMuMinus
        ptype = pp.particle.MuMinusDef()
    elseif ptype == PMuPlus
        ptype = pp.particle.MuPlusDef()
    else
        error("Type $(ptype) not supported")
    end
    propagator = pp.Propagator(ptype, joinpath(PKGDIR, "assets/proposal_config.json"))
    return propagator
end


function propagate_muon(particle; propagator=nothing, length=0)

    position = particle.position
    direction = particle.direction
    length = length > 0 ? length : particle.length
    time = particle.time
    energy = particle.energy

    propagator = isnothing(propagator) ? make_propagator(particle.type) : propagator

    initial_state = pp.particle.ParticleState()
    initial_state.energy = energy * 1E3
    initial_state.position = pp.Cartesian3D(position[1] * 100, position[2] * 100, position[3] * 100)
    initial_state.direction = pp.Cartesian3D(direction[1], direction[2], direction[3])
    initial_state.time = time / 1E9
    secondaries = propagator.propagate(initial_state, max_distance=length * 100)
    stochastic_losses = pyconvert(Vector, secondaries.stochastic_losses())
    stochastic_losses = proposal_secondary_to_particle.(stochastic_losses)

    continuous_losses = proposal_continuous_to_particle.(pyconvert(Vector, secondaries.continuous_losses()))

    T = eltype(position)

    final_state = secondaries.final_state()
    length = pyconvert(T, final_state.propagated_distance / 100)

    final_state = Particle(
        position .+ length .* direction,
        direction,
        pyconvert(T, final_state.time * 1E9),
        pyconvert(T, final_state.energy / 1E3),
        length,
        particle.type)

    return final_state, stochastic_losses, continuous_losses
end

end
