using BusInjectionModel
using Test
using HiGHS
using Ipopt
using JuMP
import OpenDSSDirect as OpenDSS


# # hack for local testing
# using Pkg
# Pkg.activate("..")
# using BusInjectionModel
# Pkg.activate(".")

CPF = BusInjectionModel.CommonOPF

SINGLE_PHASE_IEEE13_DSS_PATH = joinpath("data", "ieee13", "ieee13_makePosSeq", "Master.dss")


@testset "BusInjectionModel.jl" begin

    # @testset "IEEE13 wye only fixed point linear" begin
    #     m = Model(HiGHS.Optimizer)
    #     net = BusInjectionModel.CommonOPF.dss_to_Network(SINGLE_PHASE_IEEE13_DSS_PATH)
    #     build_bim_rectangular!(m, net, FixedPointLinear)
    # end

    @testset "IEEE13 single phase Unrelaxed model" begin
        # get the OpenDSS voltages for comparison and setting v0
        # NOTE the OpenDSS model requires lots of iterations and setting load vminpu to 0.8 to ensure
        # that the OpenDSS model converges with constant power loads
        OpenDSS.Text.Command("Redirect $SINGLE_PHASE_IEEE13_DSS_PATH")
        OpenDSS.Solution.Solve()

        @test(OpenDSS.Solution.Converged() == true)
        @test(CPF.check_opendss_powers() == true)

        dss_voltages = CPF.dss_voltages_pu()

        m = Model(Ipopt.Optimizer)
        set_optimizer_attribute(m, "print_level", 0)
        net = BusInjectionModel.CommonOPF.dss_to_Network(SINGLE_PHASE_IEEE13_DSS_PATH)
        # need to get rid of the substation regulator?
        net.v0 = dss_voltages["650"][1] + im*0
        net.Vbase = 4160 / sqrt(3)
        net.Sbase = 1e7
        net.Zbase = net.Vbase^2 / net.Sbase
        # net.bounds.v_upper_mag = net.v0 * 1.2
        # net.bounds.v_lower_mag = net.v0 * 0.7

        build_bim_rectangular!(m, net, Unrelaxed)
        M = 1e7
        @constraint(m, [t in 1:net.Ntimesteps], M >= real(m[:s0][t]) >= -M)
        @constraint(m, [t in 1:net.Ntimesteps], M >= imag(m[:s0][t]) >= -M)
        @constraint(m, [b in CPF.busses(net), t in 1:net.Ntimesteps], 
            1.1 >= imag(m[:v][b][t]) >= -1.1
        )
        @constraint(m, [b in CPF.busses(net), t in 1:net.Ntimesteps], 
            1.1 >= real(m[:v][b][t]) >= 0.8
        )
        # @constraint(m, [b in CPF.busses(net), t in 1:net.Ntimesteps], 
        #     1.2^2 >= real(m[:v][b][t])^2 + imag(m[:v][b][t])^2 >= 0.7^2
        # )
        @objective(m, Min, sum( real(m[:s0][t])^2 + imag(m[:s0][t])^2  for t in 1:net.Ntimesteps) )
        optimize!(m)
        
        @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

        r = CPF.opf_results(m, net)
        
        # remove time indices from results for convenience
        vs = Dict(b => abs.(vv[1]) for (b, vv) in pairs(r[:v]))
        
        for b in keys(vs)
            @test abs(vs[b] - dss_voltages[b][1]) < 1e-3
        end
        # NOTE bus 634 voltage is worse than others by two orders of magnitude because the
        # transformer model in CommonOPF is crude compared to the OpenDSS model.

    end

    @testset  "McCalley ISU example single phase polar voltage" begin
        # three busses in a triangle: a slack bus, a PV generator bus, and a PQ load bus

        netdict = Dict(
            :Network => Dict(:substation_bus => "1", :Sbase => 1),
            :Conductor => [
                Dict(
                    :busses => ("1", "2"),
                    :r1 => 0.0,
                    :x1 => 0.1,
                    :length => 1
                ),
                Dict(
                    :busses => ("2", "3"),
                    :r1 => 0.0,
                    :x1 => 0.1,
                    :length => 1
                ),
                Dict(
                    :busses => ("3", "1"),
                    :r1 => 0.0,
                    :x1 => 0.1,
                    :length => 1
                ),
            ],
            :Load => [
                Dict(
                    :bus => "3",
                    :kws1 => [.0028653],
                    :kvars1 => [.0012244]
                ),
            ],
            :Generator => [
                Dict(
                    :bus => "2",
                    :is_PV_bus => true,
                    :voltage_pu => [1.05],
                    :kws1 => [0.0006661]
                )
            ]
        )
        
        net = CPF.Network(netdict)
        
        m = JuMP.Model(Ipopt.Optimizer)
        set_optimizer_attribute(m, "print_level", 0)

        build_bim_polar!(m, net)

        optimize!(m)
        @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

        v_mag =  value.([(values(m[:v_mag])...)...])
        v_ang =  value.([(values(m[:v_ang])...)...])

        tol = 1e-3
        @test v_mag[2] ≈ 1.05 rtol = tol
        @test v_ang[2] * 180/pi ≈ -3 rtol = tol
        @test v_mag[3] ≈ 0.9499 rtol = tol
        @test v_ang[3] * 180/pi ≈ -10.01 rtol = tol

    end
    
end
