# PhysicsTools

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PLEnuM-group.github.io/PhysicsTools.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PLEnuM-group.github.io/PhysicsTools.jl/dev/)
[![Build Status](https://github.com/PLEnuM-group/PhysicsTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/PLEnuM-group/PhysicsTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/PLEnuM-group/PhysicsTools.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/PLEnuM-group/PhysicsTools.jl)

This package hosts various functions and types which are useful in physics-related data-science workflows:
* Interpolation: `fast_linear_interp`
* Integration: `integrate_gauss_quad`
* Coordinate transformation: `sph_to_cart`, `apply_rot`, `cart_to_sph`, `rot_to_ez_fast`, `rot_from_ez_fast`, `calc_rot_matrix`
* Probability distributions: `CategoricalSetDistribution`
* Sampling: `sample_cherenkov_track_direction`, `rand_gamma`
* Signal processing: `fwhm`
* Data manipulation: `repeat_for`, `repeat_for!`, `split_by`