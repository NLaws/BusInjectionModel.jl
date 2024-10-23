module BusInjectionModel


using CommonOPF
using JuMP
using LinearAlgebra
using Graphs

"""
    ModelType

An enum with values:
1. `FixedPointLinear`
"""
@enum ModelType begin
    FixedPointLinear  # only SinglePhase
end


export
    build_bim!,
    FixedPointLinear,
    admittance_builder,


include("single_phase_model.jl")
include("single_phase_fpl_model.jl")
include("multi_phase_fpl_model.jl")
include("utilities.jl")

end
