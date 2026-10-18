@testset "DC OPF validated against Li 2010 IEEE 5 bus" begin
    
    net = CPF.Network_IEEE5()
    m = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(m)

    build_bim_polar!(m, net, DC)

    # TODO account for cost curves
    @objective(m, Min,
        sum(
            m[:p_gen][j][gen.name][t] * gen.price_points[1] 
            for j in CPF.generator_busses(net), 
                gen in net[j][:Generator],
                t in 1:net.Ntimesteps
        )
    )

    optimize!(m)

    t = 1
    # Generator dispatch
    # Alta
    @test value(m[:p_gen]["A"]["Alta"])[t] == 40.0
    # Park City
    @test value(m[:p_gen]["A"]["Park City"])[t] == 170.0
    # Solitude
    @test isapprox(value(m[:p_gen]["C"]["Solitude"])[t], 323.49; atol = 0.01)
    # Sundance
    @test value(m[:p_gen]["D"]["Sundance"])[t] == 0.0
    # Brighton
    @test isapprox(value(m[:p_gen]["E"]["Brighton"])[t], 466.51; atol = 0.01)

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


@testset "DC OPF with battery storage" begin
    
    net = CPF.Network_IEEE5()

    net.Ntimesteps = 6
    for b in CPF.real_load_busses(net)
        load = net[b][:Load]
        load.kws1 = fill(first(load.kws1), net.Ntimesteps)
        load.kvars1 = fill(first(load.kvars1), net.Ntimesteps)
    end

    battery = CPF.Battery(
        name = "four-hour-bess",
        bus = "B",
        pmin_kw = -100.0,
        pmax_kw = 100.0,
        charge_efficiency = 0.9,
        discharge_efficiency = 0.9,
        soc_max_kwh = 400.0,
        soc_min_kwh = 0.0,
        soc_init_kwh = 0.0
    )
    
    net["B"][:Battery] = [battery]
    m = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(m)

    build_bim_polar!(m, net, DC)

    T = net.Ntimesteps

    @variables(m, begin
        z[t = 1:T], Bin
    end)


    # @constraint(m, [t in 1:T], p[t] ⟂ g[t]) 
    Mp = -1*battery.pmin_kw    # upper bound for charging
    Mg = battery.pmax_kw       # upper bound for discharging
    @constraint(m, [t in 1:T], m[:p_charge][battery.bus][battery.name][t] <= Mp * z[t])
    @constraint(m, [t in 1:T], m[:p_discharge][battery.bus][battery.name][t] <= Mg * (1 - z[t]))


    @objective(m, Min,
        sum(
            m[:p_gen][j][gen.name][t] * gen.price_points[1] 
            for j in CPF.generator_busses(net), 
                gen in net[j][:Generator],
                t in 1:net.Ntimesteps
        )
    )

    optimize!(m)

    @test JuMP.termination_status(m) == JuMP.OPTIMAL
    # no BESS usage with flat offers and load
    @test isapprox(
        JuMP.value.(m[:p_charge]["B"][battery.name]),
        zeros(net.Ntimesteps);
        atol = 1e-7,
    )
    @test isapprox(
        JuMP.value.(m[:p_discharge]["B"][battery.name]),
        zeros(net.Ntimesteps);
        atol = 1e-7,
    )

    # now change a generators offer curve to incentivize BESS usage
    # do we have time varying gen offers yet?
end


@testset "DC Copper Plate" begin
    net = BusInjectionModel.CommonOPF.Network(COPPER_PLATE_MODEL_PATH)
    m = JuMP.Model(HiGHS.Optimizer)
    JuMP.set_silent(m)
    build_bim_polar!(m, net, DC)
    thermal_price = 50.0
    renewable_price = 0.0
    charging_cost = 1.0
    discharging_cost = 1.0

    # perfect foresight bids and offers
    ess_bids = [43.1, 50.0, 0.0, 0.0, 50.0, 40.0]
    ess_offers = [43.1, 50.0, 0.0, 0.0, 49.9, 49.9]

    # TODO getter for renewable generators and setting time varying limits
    @constraint(m, [t in 1:net.Ntimesteps], m[:p_gen]["b1"]["renewable"][t] <= net["b1"][:Generator][1].kws1[t])

    # TODO CommonOPF.Generator.pmax -> pmax_kw, pmin -> pmin_kw
    @objective(m, Min,
        (1/net.time_steps_per_hour) * (
            sum(m[:p_gen]["b1"]["renewable"][t] for t in 1:net.Ntimesteps) * renewable_price +
            sum(m[:p_gen]["b1"]["thermal"][t] for t in 1:net.Ntimesteps)   * thermal_price   +
            # degradation cost
            sum(m[:p_charge]["b1"]["battery"][t] for t in 1:net.Ntimesteps) * charging_cost +
            sum(m[:p_discharge]["b1"]["battery"][t] for t in 1:net.Ntimesteps) * discharging_cost
        ) - 
        # battery energy value
        40 * (m[:soc]["b1"]["battery"][net.Ntimesteps] - m[:soc]["b1"]["battery"][0])
        # m[:p_discharge]["b1"]["battery"][t] * ess_offers[t] -
        # m[:p_charge]["b1"]["battery"][t]    * ess_bids[t]
    )

    optimize!(m)

    @test JuMP.termination_status(m) == JuMP.OPTIMAL

    prices = value.(dual.((m[:load_balance]))) * net.time_steps_per_hour
    discharge = value.(m[:p_discharge]["b1"]["battery"])
    charge = value.(m[:p_charge]["b1"]["battery"])
    thermal = value.(m[:p_gen]["b1"]["thermal"])
    renewable = value.(m[:p_gen]["b1"]["renewable"])
    soc = value.(m[:soc]["b1"]["battery"])
    terminal_soc_value = 40 * (soc[net.Ntimesteps] - soc[0])

    # The load-balance dual is in cost per model power unit. Convert it to an
    # energy price using the same time scaling as the objective.
    @test isapprox(prices, [38.69, 50.0, -0.0, -0.0, 50.0, 50.0], atol = 1e-2)

    ess_surplus = discharge' * prices - charge' * prices -
        # degradation cost
        charging_cost * sum(charge) - discharging_cost * sum(discharge)
    ess_surplus = round(ess_surplus / net.time_steps_per_hour + terminal_soc_value, digits = 2)
    # Unlike print_results in the reference implementation, this calculation
    # uses unrounded solver values, so the expected surplus rounds to 34391.11.
    @test isapprox(ess_surplus, 34391.11, atol = 1e-2)

    # Charging and renewable dispatch are unique for this case.
    @test isapprox(charge, [100.0, 0.0, 200.0, 22.2222, 0.0, 0.0], atol = 1e-2)
    @test isapprox(renewable, [200.0, 100.0, 500.0, 322.2222, 25.0, 0.0], atol = 1e-2)

    # Discharge and thermal dispatch are unique through period 4.
    @test isapprox(discharge[1:4], [0.0, 81.0, 0.0, 0.0], atol = 1e-2)
    @test isapprox(thermal[1:4], [0.0, 19.0, 0.0, 0.0], atol = 1e-2)

    # Periods 5 and 6 have the same energy price and time-invariant operating
    # costs. Shifting discharge and thermal generation between those periods
    # leaves the objective unchanged, so their individual dispatch is not
    # unique. Test the invariant totals instead of a solver-specific split.
    @test isapprox(sum(discharge[5:6]), 180.0, atol = 1e-2)
    @test isapprox(sum(thermal[5:6]), 95.0, atol = 1e-2)
    @test isapprox(sum(discharge), 261.0, atol = 1e-2)
    @test isapprox(sum(thermal), 114.0, atol = 1e-2)

    # The period-5 SOC depends on the non-unique dispatch split above. SOC is
    # unique through period 4, and every optimum fully discharges by period 6.
    @test isapprox([soc[t] for t in 1:4], [360.0, 0.0, 720.0, 800.0], atol = 1e-2)
    @test isapprox(soc[net.Ntimesteps], 0.0, atol = 1e-2)
    @test isapprox(objective_value(m), 25132.89, atol = 1e-2)
    
end
