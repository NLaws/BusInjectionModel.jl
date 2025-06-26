@testset "multiphase fixed point" begin

    m = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(m)
    net = CPF.Network_IEEE13()

    build_bim_rectangular!(m, net, FixedPointLinear)

    # using Logging
    # global_logger(ConsoleLogger(stderr, Logging.Debug))

    BusInjectionModel.solve_fixed_point_to_tol(m, net)

    v_fp = abs.(value.(m[:v]).data)
    v_fp_by_bus = Dict()
    for (term, val) in zip(m[:ll_terminals], v_fp)
        if !(term.bus in keys(v_fp_by_bus))
            v_fp_by_bus[term.bus] = [0.0, 0.0, 0.0]
        end
        v_fp_by_bus[term.bus][term.phase] = val
    end

    # compare to the unrelaxed solution
    m2 = Model(Ipopt.Optimizer)
    JuMP.set_silent(m)

    build_bim_rectangular!(m2, net, Unrelaxed)

    M = 1e3
    @constraint(m2, [t in 1:net.Ntimesteps], M .>= real(m2[:s0][net.substation_bus][t]) .>= -M)

    @constraint(m2, [t in 1:net.Ntimesteps], M .>= imag(m2[:s0][net.substation_bus][t]) .>= -M)

    non_sub_busses = setdiff(CPF.busses(net), [net.substation_bus])
    # NOTE that the real and imaginary voltage parts can be negative
    @constraint(m2, 
        [b in non_sub_busses, t in 1:net.Ntimesteps, phs in CPF.phases_connected_to_bus(net, b)], 
        1.1 >= imag(m2[:v][b][t][phs]) >= -1.1
    )
    @constraint(m2, 
        [b in non_sub_busses, t in 1:net.Ntimesteps, phs in CPF.phases_connected_to_bus(net, b)],
        1.5 .>= real(m2[:v][b][t][phs]) .>= -1.1
    )

    @objective(m2, Min, sum( 
        real(m2[:s0][net.substation_bus][t][phs]).^2 
        + imag(m2[:s0][net.substation_bus][t][phs]).^2  
        for t in 1:net.Ntimesteps, phs in 1:3) 
    )

    optimize!(m2)
    r = CPF.opf_results(m2, net)
        
    # remove time indices from results for convenience
    vs = Dict(b => abs.(vv[1]) for (b, vv) in pairs(r[:v]))
        
    for (bus, vals) in vs
        for (phs, val) in enumerate(vals)
            # no substation_bus in v_fp_by_bus
            if !(bus == net.substation_bus)
                @test v_fp_by_bus[bus][phs] ≈ val rtol=1e-3
            end
        end
    end


end
