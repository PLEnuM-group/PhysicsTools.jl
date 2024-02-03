using PhysicsTools
using Test
using DataStructures
using StaticArrays
using Distributions

@testset "PhysicsTools.jl" begin
    @testset "Utils" begin

        @testset "repeat_for" begin
            let x = [[1, 2, 3] [4, 5, 6] [7, 8, 9]], y = [[0, 0, 0] [0, 0, 0] [0, 0, 0] [0, 0, 0] [0, 0, 0] [0, 0, 0]]
                @test repeat_for(x, [1, 2, 3]) == [[1, 2, 3] [4, 5, 6] [4, 5, 6] [7, 8, 9] [7, 8, 9] [7, 8, 9]]
                @test repeat_for(x, [0, 2, 0]) == [[4, 5, 6] [4, 5, 6]]
                @test repeat_for(x, [0, 0, 0]) == reshape([[] [] []], 3, 0)
                @test repeat_for!(x, [1, 2, 3], y) == [[1, 2, 3] [4, 5, 6] [4, 5, 6] [7, 8, 9] [7, 8, 9] [7, 8, 9]]
            end
            let x = [["a", "b", "c"] ["d", "e", "f"] ["g", "h", "i"]]
                @test repeat_for(x, [1, 2, 3]) == [["a", "b", "c"] ["d", "e", "f"] ["d", "e", "f"] ["g", "h", "i"] ["g", "h", "i"] ["g", "h", "i"]]
            end
        end

        @testset "split_by" begin
            let x = [1, 2, 3, 4, 5, 6, 7, 8, 9], y = [[], []]
                @test split_by(x, [2, 3, 4]) == [[1, 2], [3, 4, 5], [6, 7, 8, 9]]
                @test split_by(x, [2, 0, 1]) == [[1, 2], [], [3]]
                @test split_by(x, [1, 2]) == [[1], [2, 3]]
                @test split_by(x, [0, 0]) == [[], []]
                @test split_by!(x, [1, 2], y) == [[1], [2, 3]]
                #handling out of bounds indicies
            end
            let x = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
                @test split_by(x, [1, 2]) == [["a"], ["b", "c"]]
            end
        end

        @testset "sph_to_cart" begin
            let theta = 0.0, phi = 0.0
                @test sph_to_cart([theta, phi]) ≈ SA[0.0, 0.0, 1.0]
            end

            let theta = 0.0, phi = π
                @test sph_to_cart([theta, phi]) ≈ SA[0.0, 0.0, 1.0]
            end

            let theta = π / 2, phi = 0
                @test sph_to_cart([theta, phi]) ≈ SA[1.0, 0.0, 0.0]
            end

            let theta = π / 2, phi = π
                @test sph_to_cart([theta, phi]) ≈ SA[-1.0, 0.0, 0.0]
            end

            let theta = π / 2, phi = π / 2
                @test sph_to_cart([theta, phi]) ≈ SA[0.0, 1.0, 0.0]
            end
        end

        @testset "cyl_to_cart" begin
            let theta = 1, phi = 0.0, z = 1
                @test cyl_to_cart(theta, phi, z) ≈ SA[1.0, 0.0, 1.0]
            end
            let theta = 1, phi = π, z = 1
                @test cyl_to_cart(theta, phi, z) ≈ SA[-1.0, 0.0, 1.0]
            end
            let theta = 1, phi = π / 2, z = 1
                @test cyl_to_cart(theta, phi, z) ≈ SA[0.0, 1.0, 1.0]
            end
            let theta = 1, phi = 3 * π / 2, z = -1
                @test cyl_to_cart(theta, phi, z) ≈ SA[0.0, -1.0, -1.0]
            end
        end

        @testset "cyl_to_cartV" begin
            let theta = 1, phi = 0.0, z = 1
                @test cyl_to_cart([theta, phi, z]) ≈ SA[1.0, 0.0, 1.0]
            end
            let theta = 1, phi = π, z = 1
                @test cyl_to_cart([theta, phi, z]) ≈ SA[-1.0, 0.0, 1.0]
            end
            let theta = 1, phi = π / 2, z = 1
                @test cyl_to_cart([theta, phi, z]) ≈ SA[0.0, 1.0, 1.0]
            end
            let theta = 1, phi = 3 * π / 2, z = -1
                @test cyl_to_cart([theta, phi, z]) ≈ SA[0.0, -1.0, -1.0]
            end
        end

        @testset "cart_to_cyl" begin
            let theta = 0.3, phi = 1.5, z = 1
                @test all(cart_to_cyl(cyl_to_cart(theta, phi, z)) .≈ (theta, phi, z))
            end
            let theta = π, phi = 1, z = 1
                @test all(cart_to_cyl(cyl_to_cart(theta, phi, z)) .≈ (theta, phi, z))
            end
            let theta = 0.3, phi = π / 2, z = -1
                @test all(cart_to_cyl(cyl_to_cart(theta, phi, z)) .≈ (theta, phi, z))
            end
            let theta = 2.5, phi = 3.2, z = -1
                @test all(cart_to_cyl(cyl_to_cart(theta, phi, z)) .≈ (theta, phi, z))
            end
            let theta = 1.5, phi = 0.1, z = -1
                @test all(cart_to_cyl(cyl_to_cart(theta, phi, z)) .≈ (theta, phi, z))
            end
        end

        #@testset "calc_gamma_shape_mean_fwhm" begin
        #    let mean=4, tfwhm=1
        #        @test calc_gamma_shape_mean_fwhm(4, 1)
        #    end
        #end

    end
end


#How to test errors?
#chech negative numbers
#check out of bound indicies
#cart to cyl imprecise 0,0,1