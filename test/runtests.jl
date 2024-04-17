using PhysicsTools
using Test
using DataStructures
using StaticArrays
using Distributions

@testset "PhysicsTools.jl" begin

    include("runtestsParticleA.jl")
    include("runtestsUtilsA.jl")

    @testset "Utils" begin

        @testset "fast_linear_interp" begin
            # test edge cases
            let knots = [1, 2, 3], lower = 1, upper = 3
                x_eval = 1
                @test fast_linear_interp(x_eval, knots, lower, upper) == knots[1]
                x_eval = 3
                @test fast_linear_interp(x_eval, knots, lower, upper) == knots[3]
                x_eval = 0
                @test fast_linear_interp(x_eval, knots, lower, upper) == knots[1]
                x_eval = 4
                @test fast_linear_interp(x_eval, knots, lower, upper) == knots[3]
            end

            let xs = [1, 2, 3], ys = [10, 20, 30]
                x_eval = 1
                @test fast_linear_interp(x_eval, xs, ys) == ys[1]
                x_eval = 3
                @test fast_linear_interp(x_eval, xs, ys) == ys[3]
                x_eval = 0
                @test fast_linear_interp(x_eval, xs, ys) == ys[1]
                x_eval = 4
                @test fast_linear_interp(x_eval, xs, ys) == ys[3]
            end

            # test normal cases
            let knots = [1, 5, 7], lower = 1, upper = 3
                x_eval = 2.5
                @test fast_linear_interp(x_eval, knots, lower, upper) ≈ 6
            end

            let xs = [1, 2, 3], ys = [10, 20, 30]
                x_eval = 2.5
                @test fast_linear_interp(x_eval, xs, ys) ≈ 25
            end

            # test type stability

            let knots = [1, 5, 7], lower = 1, upper = 3
                x_eval = 2.5
                @inferred fast_linear_interp(x_eval, knots, lower, upper)
            end

            let knots = [1, 5, 7], lower = 1, upper = 3
                x_eval = 2
                @inferred fast_linear_interp(x_eval, knots, lower, upper)
            end

            let knots = [1.0, 5.0, 7.0], lower = 1.0, upper = 3.0
                x_eval = 2
                @inferred fast_linear_interp(x_eval, knots, lower, upper)
            end

            # test Float32
            let knots = [1.0f0, 5.0f0, 7.0f0], lower = 1.0f0, upper = 3.0f0
                x_eval = 2
                typeof(fast_linear_interp(x_eval, knots, lower, upper)) == Float32
            end


        end


        @testset "Rotations" begin

            e_z = SA[0.0, 0.0, 1.0]
            e_x = SA[1.0, 0.0, 0.0]
            e_y = SA[0.0, 1.0, 0.0]
            # Test that applying the rotation to the first vector yields the second vector
            let theta = 0.2, phi = 1.4
                dir = sph_to_cart(theta, phi)
                @test all(isapprox.(apply_rot(dir, e_z, dir), e_z; atol=1E-9))
            end

            let θ_1 = π / 2, θ_2 = π / 2, ϕ_1 = 0, ϕ_2 = π / 4
                # Rotation by 45° around z
                dir_1 = sph_to_cart(θ_1, ϕ_1)
                dir_2 = sph_to_cart(θ_2, ϕ_2)
                @test all(isapprox.(apply_rot(dir_1, dir_2, e_x), SA[1/sqrt(2), 1/sqrt(2), 0]; atol=1E-9))
            end

            let θ_1 = 0.2, θ_2 = 1.3, ϕ_1 = 1.8, ϕ_2 = 2.5,
                dir_1 = sph_to_cart(θ_1, ϕ_1)

                dir_2 = sph_to_cart(θ_2, ϕ_2)
                @test all(apply_rot(dir_1, e_z, dir_2) .≈ rot_to_ez_fast(dir_1, dir_2))
                @test all(apply_rot(e_z, dir_1, dir_2) .≈ rot_from_ez_fast(dir_1, dir_2))

            end

        end

        @testset "integrate_gauss_quad" begin
            let f = x -> x^3
                @test integrate_gauss_quad(f, 0.0, 1.0) ≈ 1 / 4
            end

            let f = cos, a = 0.0, b = 0.5
                @test integrate_gauss_quad(f, a, b) ≈ sin(b) - sin(a)
            end

        end

        @testset "sph_to_cart" begin
            let theta = 0.0, phi = 0.0
                @test sph_to_cart(theta, phi) ≈ SA[0.0, 0.0, 1.0]
            end

            let theta = 0.0, phi = π
                @test sph_to_cart(theta, phi) ≈ SA[0.0, 0.0, 1.0]
            end

            let theta = π / 2, phi = 0
                @test sph_to_cart(theta, phi) ≈ SA[1.0, 0.0, 0.0]
            end

            let theta = π / 2, phi = π
                @test sph_to_cart(theta, phi) ≈ SA[-1.0, 0.0, 0.0]
            end

            let theta = π / 2, phi = π / 2
                @test sph_to_cart(theta, phi) ≈ SA[0.0, 1.0, 0.0]
            end
        end

        @testset "cart_to_sph" begin

            let theta = 0.3, phi = 1.5
                @test all(cart_to_sph(sph_to_cart(theta, phi)) .≈ (theta, phi))
            end

            let theta = 0, phi = 1
                @test all(cart_to_sph(sph_to_cart(theta, phi)) .≈ (0, 0))
            end

            let theta = π, phi = 1
                @test all(cart_to_sph(sph_to_cart(theta, phi)) .≈ (π, 0))
            end


            let theta = 0.3, phi = π / 2
                @test all(cart_to_sph(sph_to_cart(theta, phi)) .≈ (theta, phi))
            end

            let theta = 2.5, phi = 3.2
                @test all(cart_to_sph(sph_to_cart(theta, phi)) .≈ (theta, phi))
            end

            let theta = 1.5, phi = 0.1
                @test all(cart_to_sph(sph_to_cart(theta, phi)) .≈ (theta, phi))
            end

        end

        @testset "CategoricalSetDistribution" begin
            let pdist = CategoricalSetDistribution(OrderedSet([:EMinus, :EPlus]), [1.0, 0.0])
                @test rand(pdist) === :EMinus
            end

            let pdist = CategoricalSetDistribution(OrderedSet([:EMinus, :EPlus]), Categorical([1.0, 0.0]))
                @test rand(pdist) === :EMinus
            end

            let err = nothing
                try
                    CategoricalSetDistribution(OrderedSet([:EMinus, :EPlus]), Categorical([1,]))
                catch err
                end
                @test err isa Exception
            end

            let err = nothing
                try
                    CategoricalSetDistribution(OrderedSet[:EMinus, :EPlus], [1.0])
                catch err
                end
                @test err isa Exception
            end
        end
    end
end
