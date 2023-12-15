module ProposalInterface
using PyCall
using StaticArrays
using Conda
using Logging
import Pkg

import ..Particle, ..PEMinus, ..PEPlus, ..PMuMinus, ..PMuPlus, ..PHadronShower, ..ParticleType


export proposal_secondary_to_particle, propagate_muon

const pp = PyNULL()

const PKGDIR = pkgdir(@__MODULE__)

function __init__()
    
    try
        tmp = pyimport("proposal")
    catch
        @warn "Could not import proposal. Attempting to install..."
        ENV["PYTHON"] = ""
        Pkg.build("PyCall")
        Conda.pip_interop(true)
        Conda.pip("install", "proposal")
        tmp = pyimport("proposal")
    end
   
    if !isdir(joinpath(PKGDIR, "assets/proposal_tables"))
        mkpath(joinpath(PKGDIR, "assets/proposal_tables"))
    end
    tmp.InterpolationSettings.tables_path = joinpath(PKGDIR, "assets/proposal_tables")
    copy!(pp, tmp)
end

function proposal_secondary_to_particle(loss)
    energy = loss.energy / 1E3
    pos = SA[loss.position.x/100, loss.position.y/100, loss.position.z/100]
    dir = SA[loss.direction.x, loss.direction.y, loss.direction.z]
    time = loss.time * 1E9

    pname = pp.particle.Interaction_Type(loss.type).name

    if pname == "photonuclear"
        ptype = PHadronShower
    else
        ptype = PEMinus
    end

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


function propagate_muon(particle)

    position = particle.position
    direction = particle.direction
    length = particle.length
    time = particle.time
    energy = particle.energy

    propagator = make_propagator(particle.type)

    initial_state = pp.particle.ParticleState()
    initial_state.energy = energy * 1E3
    initial_state.position = pp.Cartesian3D(position[1] * 100, position[2] * 100, position[3] * 100)
    initial_state.direction = pp.Cartesian3D(direction[1], direction[2], direction[3])
    initial_state.time = time / 1E9
    secondaries = propagator.propagate(initial_state, max_distance=length * 100)
    stochastic_losses = secondaries.stochastic_losses()
    stochastic_losses = proposal_secondary_to_particle.(stochastic_losses)


    T = eltype(position)

    final_state = secondaries.final_state()
    length = T(final_state.propagated_distance / 100)

    final_state = Particle(
        position .+ length .* direction,
        direction,
        T(final_state.time * 1E9),
        T(final_state.energy / 1E3),
        length,
        particle.type)

    return final_state, stochastic_losses
end

end
