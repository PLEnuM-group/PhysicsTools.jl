module Types
using StaticArrays

include("particle.jl")

export Position, Direction

const Position{T} = SVector{3, T}
const Direction{T} = SVector{3, T}

end