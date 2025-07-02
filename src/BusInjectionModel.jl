module BusInjectionModel


using CommonOPF
using JuMP
using LinearAlgebra
using Printf


"""
    ModelType

An enum with values:
1. `FixedPointLinear`
2. `Unrelaxed`
"""
@enum ModelType begin
    FixedPointLinear
    Unrelaxed
end


export
    build_bim_polar!,
    build_bim_rectangular!,
    solve_fixed_point_to_tol,
    FixedPointLinear,
    Unrelaxed


include("single_phase_model.jl")
include("multi_phase_model.jl")
include("fixed_point_linear_model.jl")

end
