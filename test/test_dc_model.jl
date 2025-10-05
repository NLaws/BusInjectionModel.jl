@testset "DC OPF validated against Li 2019 IEEE 5 bus" begin
    
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

    # LMPs TODO

    # Line Flows TODO

end