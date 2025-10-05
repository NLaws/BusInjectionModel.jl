@testset "DC OPF validated against Li 2010 IEEE 5 bus" begin
    
    net = CPF.Network_IEEE5()
    m = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(m)

    build_bim_polar!(m, net, DC)

    # TODO account for cost curves
    @objective(m, Min,
        sum(
            m[:p_gen][j][i][t] * gen.price_points[1] 
            for j in CPF.generator_busses(net), 
                (i,gen) in enumerate(net[j][:Generator]), 
                t in 1:net.Ntimesteps
        )
    )

    optimize!(m)

    t = 1
    # Generator dispatch
    # Alta
    @test value(m[:p_gen]["A"][1])[t] == 40.0
    # Park City
    @test value(m[:p_gen]["A"][2])[t] == 170.0
    # Solitude
    @test isapprox(value(m[:p_gen]["C"][1])[t], 323.49; atol = 0.01)
    # Sundance
    @test value(m[:p_gen]["D"][1])[t] == 0.0
    # Brighton
    @test isapprox(value(m[:p_gen]["E"][1])[t], 466.51; atol = 0.01)

    # LMPs
    congestion_costs = JuMP.dual.(m[:bus_real_power_injection_constraints]...)
    LMP_energy = JuMP.dual(m[:load_balance][t])
    LMPs = Dict()
    for b in CPF.busses(net)
        LMPs[b] = JuMP.dual(m[:pj_constraints][b][1])
    end

    @test isapprox(LMPs["A"], 16.98; atol = 0.01)
    @test isapprox(LMPs["B"], 26.38; atol = 0.01)
    @test isapprox(LMPs["C"], 30.00; atol = 0.01)
    @test isapprox(LMPs["D"], 39.94; atol = 0.01)
    @test isapprox(LMPs["E"], 10.00; atol = 0.01)

    @test isapprox(congestion_costs["A"], -22.97; atol = 0.01)
    @test isapprox(congestion_costs["B"], -13.56; atol = 0.01)
    @test isapprox(congestion_costs["C"], -9.94; atol = 0.01)
    @test isapprox(congestion_costs["D"],  0.00; atol = 0.01)
    @test isapprox(congestion_costs["E"], -29.94; atol = 0.01)

    @test isapprox(LMP_energy, 39.94; atol=0.01)

    # Line Flows
    line_flows = Dict()
    for e in CPF.edges(net)
        line_flows[e] = value(m[:line_flow][e][t])
    end

    # NOTE signs reflect the node order in BusInjectionModel, not Li 2019
    # Flow directions have been confirmed manually
    @test isapprox(line_flows[("B", "A")], -249.72; atol = 0.01)
    @test isapprox(line_flows[("A", "D")],  186.79; atol = 0.01)
    @test isapprox(line_flows[("A", "E")], -226.51; atol = 0.01)
    @test isapprox(line_flows[("B", "C")], -50.28; atol = 0.01)
    @test isapprox(line_flows[("C", "D")], -26.79; atol = 0.01)
    @test isapprox(line_flows[("D", "E")], -240.00; atol = 0.01)

end