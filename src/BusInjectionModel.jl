module BusInjectionModel


using CommonOPF
using JuMP
using LinearAlgebra


"""
    ModelType

An enum with values:
1. `FixedPointLinear`
"""
@enum ModelType begin
    FixedPointLinear  # only SinglePhase
    Unrelaxed
end


export
    build_bim!,
    FixedPointLinear,
    Unrelaxed


include("single_phase_model.jl")
include("single_phase_fpl_model.jl")

end
