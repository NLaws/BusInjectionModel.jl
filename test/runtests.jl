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
MULTIPHASE_IEEE13_DSS_PATH = joinpath("data", "ieee13", "IEEE13Nodeckt_no_trfxs.dss")


@testset "BusInjectionModel.jl" begin

    include("test_fixed_point.jl")

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

    @testset  "McCalley ISU example: single-phase, 3-bus, polar voltage" begin
        # three busses in a triangle: a slack bus, a PV generator bus, and a PQ load bus

        netdict = Dict(
            :Network => Dict(
                :substation_bus => "1", 
                :Sbase => 1,
                :v_upper_mag => 1.1,
                :v_lower_mag => 0.9,
                :v_upper_ang => π,
                :v_lower_ang => -π,
            ),
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
                    :voltage_series_pu => [1.05],
                    :kws1 => [0.0006661]
                )
            ]
        )
        
        net = CPF.Network(netdict)
        # TODO these four tests should be in CommonOPF
        @test net.bounds.v_upper_mag == 1.1
        @test net.bounds.v_lower_mag == 0.9
        @test net.bounds.v_upper_ang == π
        @test net.bounds.v_lower_ang == -π
        
        m = JuMP.Model(Ipopt.Optimizer)
        set_optimizer_attribute(m, "print_level", 0)

        build_bim_polar!(m, net)

        optimize!(m)
        @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

        v_mag =  value.([(values(m[:v_mag])...)...])
        v_ang =  value.([(values(m[:v_ang])...)...])
        p =  value.([(values(m[:pj])...)...])
        q =  value.([(values(m[:qj])...)...])

        tol = 1e-3
        @test v_mag[2] ≈ 1.05 rtol = tol
        @test v_ang[2] * 180/pi ≈ -3 rtol = tol
        @test v_mag[3] ≈ 0.9499 rtol = tol
        @test v_ang[3] * 180/pi ≈ -10.01 rtol = tol

        @test p[2] ≈ net["2"][:Generator].kws1[1] * 1e3
        @test q[3] ≈ -net["3"][:Load].kvars1[1] * 1e3

    end
    
    @testset "IEEE13 multiphase Unrelaxed rectangular model" begin

        # get the OpenDSS voltages for comparison and setting v0
        OpenDSS.Text.Command("Redirect $MULTIPHASE_IEEE13_DSS_PATH")
        OpenDSS.Solution.Solve()

        @test(OpenDSS.Solution.Converged() == true)
        @test(CPF.check_opendss_powers() == true)

        dss_voltages = CPF.dss_voltages_pu()

        net = CPF.dss_to_Network(MULTIPHASE_IEEE13_DSS_PATH)
        net.v0 = dss_voltages["650"][1]
        net.Vbase = 4160 / sqrt(3)
        net.Sbase = 1e6
        net.Zbase = net.Vbase^2 / net.Sbase

        m = Model(Ipopt.Optimizer)
        set_optimizer_attribute(m, "print_level", 0)


        build_bim_rectangular!(m, net, Unrelaxed)

        # add bounds
        # TODO more automated bounds and good defaults
        M = 1e3
        @constraint(m, [t in 1:net.Ntimesteps], M .>= real(m[:s0][net.substation_bus][t]) .>= -M)

        @constraint(m, [t in 1:net.Ntimesteps], M .>= imag(m[:s0][net.substation_bus][t]) .>= -M)

        non_sub_busses = setdiff(CPF.busses(net), [net.substation_bus])
        # NOTE that the real and imaginary voltage parts can be negative
        @constraint(m, 
            [b in non_sub_busses, t in 1:net.Ntimesteps, phs in CPF.phases_connected_to_bus(net, b)], 
            1.1 >= imag(m[:v][b][t][phs]) >= -1.1
        )
        @constraint(m, 
            [b in non_sub_busses, t in 1:net.Ntimesteps, phs in CPF.phases_connected_to_bus(net, b)],
            1.5 .>= real(m[:v][b][t][phs]) .>= -1.1
        )

        @objective(m, Min, sum( 
            real(m[:s0][net.substation_bus][t][phs]).^2 
            + imag(m[:s0][net.substation_bus][t][phs]).^2  
            for t in 1:net.Ntimesteps, phs in 1:3) 
        )

        optimize!(m)
        @test termination_status(m) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED]

        ## for inspecting values
        # vs = value.([(((values(values(m[:v]))...)...)...)...])
        # v_mags = abs.(vs)
        # s0 = value.(m[:s0][net.substation_bus][1])

        r = CPF.opf_results(m, net)
        
        # remove time indices from results for convenience
        vs = Dict(b => abs.(vv[1]) for (b, vv) in pairs(r[:v]))
        
        for b in keys(vs)
            for (i, phsv) in enumerate(filter(v -> v != 0, vs[b]))
                @test abs(phsv - dss_voltages[b][i]) < 2e-4
            end
        end
        
    end

    # @testset "IEEE 118 Single Phase" begin
    # TODO why adding gen real power constraint makes IEEE 118 problem infeasible?
    # can include gen real power if remove the gen v_mag constraints, but then get nutty voltages
    #     net = CPF.Network_IEEE118()
    #     # m = JuMP.Model(HiGHS.Optimizer)
    #     m = JuMP.Model(Ipopt.Optimizer)
    #     # JuMP.set_silent(m)
    #     # build_bim_rectangular!(m, net, FixedPointLinear)

    #     net.bounds.v_upper_mag = 1.1
    #     net.bounds.v_lower_mag = 0.8
    #     net.bounds.v_upper_ang = π
    #     net.bounds.v_lower_ang = -π
    #     build_bim_polar!(m, net, Unrelaxed)

    #     optimize!(m)

    #     v_mag =  value.([(values(m[:v_mag])...)...])

    #     pj =  value.([(values(m[:pj])...)...])

    #     q_gen =  value.([(values(m[:q_gen])...)...])

    #     # @objective(m, Min, sum(real(m[:v][term, t]) for term in m[:ll_terminals], t in 1:net.Ntimesteps))

    #     # using Logging
    #     # global_logger(ConsoleLogger(stderr, Logging.Debug))

    #     # BusInjectionModel.solve_fixed_point_to_tol(m, net)

    # end

#     @testset "IEEE 8500 Node" begin
#         # PROBLEM IS INFEASIBLE
#         net = CPF.Network_IEEE8500()
#         net.Vbase = 12470 / sqrt(3)
#         net.Sbase = 1e9
#         net.Zbase = net.Vbase^2 / net.Sbase


#         if !(@isdefined GRB_ENV)  # the whole point is to only create GRB_ENV once per Julia session
#             const GRB_ENV = Gurobi.Env()
#         end
#         m = Model(() -> Gurobi.Optimizer(GRB_ENV))

#         m = Model(Ipopt.Optimizer)

#         build_bim_rectangular!(m, net, Unrelaxed)


#         # add bounds
#         # TODO more automated bounds and good defaults
#         M = 1e9
#         @constraint(m, [t in 1:net.Ntimesteps], M .>= real(m[:s0][net.substation_bus][t]) .>= -M)

#         @constraint(m, [t in 1:net.Ntimesteps], M .>= imag(m[:s0][net.substation_bus][t]) .>= -M)

#         non_sub_busses = setdiff(CPF.busses(net), [net.substation_bus])
#         # NOTE that the real and imaginary voltage parts can be negative
#         @constraint(m, 
#             [b in non_sub_busses, t in 1:net.Ntimesteps, phs in CPF.phases_connected_to_bus(net, b)], 
#             1.1 >= imag(m[:v][b][t][phs]) >= -1.1
#         )
#         @constraint(m, 
#             [b in non_sub_busses, t in 1:net.Ntimesteps, phs in CPF.phases_connected_to_bus(net, b)],
#             1.5 .>= real(m[:v][b][t][phs]) .>= -1.1
#         )

#         @objective(m, Min, sum( 
#             real(m[:s0][net.substation_bus][t][phs]).^2 
#             + imag(m[:s0][net.substation_bus][t][phs]).^2  
#             for t in 1:net.Ntimesteps, phs in 1:3) 
#         )

#         optimize!(m)

#     end
end