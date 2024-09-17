using BusInjectionModel
using Test
using HiGHS
using JuMP


# # hack for local testing
# using Pkg
# Pkg.activate("..")
# using LinearPowerFlow
# Pkg.activate(".")


@testset "BusInjectionModel.jl" begin

    @testset "IEEE13 wye only" begin
        m = Model(HiGHS.Optimizer)
        dss_path = joinpath("data", "ieee13", "IEEE13Nodeckt_no_trfxs.dss")
        net = BusInjectionModel.CommonOPF.dss_to_Network(dss_path)
        # build_bim!(m, net, FixedPointLinear)
    end
    
end
