using PhysicsTools
using Test
using DataStructures
using StaticArrays
using Distributions

@testset "PhysicsTools.jl" begin
    @testset "particle" begin

        let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PEMinus)
            @test p.position == SVector(1, 1, 1)
            @test p.direction == SVector(2, 2, 2)
            @test p.time == 3
            @test p.energy == 4
            @test p.length == 5
            @test p.type == PEMinus
        end

        @testset "pdg_codes" begin
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PEPlus)
                @test pdg_code(p.type) == -11
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PEMinus)
                @test pdg_code(p.type) == 11
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PGamma)
                @test pdg_code(p.type) == 22
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PMuMinus)
                @test pdg_code(p.type) == 13
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PMuPlus)
                @test pdg_code(p.type) == -13
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuE)
                @test pdg_code(p.type) == 12
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuEBar)
                @test pdg_code(p.type) == -12
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuMu)
                @test pdg_code(p.type) == 14
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuMuBar)
                @test pdg_code(p.type) == -14
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuTau)
                @test pdg_code(p.type) == 16
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuTauBar)
                @test pdg_code(p.type) == -16
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PUnknown)
                @test pdg_code(p.type) == 0
            end
        end

        @testset "ptype_for_code" begin
            @test ptype_for_code(-11) == PEPlus
            @test ptype_for_code(11) == PEMinus
            @test ptype_for_code(22) == PGamma
            @test ptype_for_code(13) == PMuMinus
            @test ptype_for_code(-13) == PMuPlus
            @test ptype_for_code(12) == PNuE
            @test ptype_for_code(-12) == PNuEBar
            @test ptype_for_code(14) == PNuMu
            @test ptype_for_code(-14) == PNuMuBar
            @test ptype_for_code(16) == PNuTau
            @test ptype_for_code(-16) == PNuTauBar
            @test ptype_for_code(0) == PUnknown
        end

        @testset "is_neutrino" begin
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PEPlus)
                @test is_neutrino(p) == false
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuE)
                @test is_neutrino(p) == true
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuEBar)
                @test is_neutrino(p) == true
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuMu)
                @test is_neutrino(p) == true
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuMuBar)
                @test is_neutrino(p) == true
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuTau)
                @test is_neutrino(p) == true
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PNuTauBar)
                @test is_neutrino(p) == true
            end
        end

        @testset "particle_shape" begin
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PEPlus)
                @test particle_shape(p) == Cascade()
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PEMinus)
                @test particle_shape(p) == Cascade()
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PGamma)
                @test particle_shape(p) == Cascade()
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PMuMinus)
                @test particle_shape(p) == Track()
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PMuPlus)
                @test particle_shape(p) == Track()
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PHadronShower)
                @test particle_shape(p) == Cascade()
            end
            let p = Particle(SVector(1, 1, 1), SVector(2, 2, 2), 3, 4, 5, PLightSabre)
                @test particle_shape(p) == Track()
            end
        end

    end
end

#No PUnknown in ALL_PARTICLES
#How to test Show() functions