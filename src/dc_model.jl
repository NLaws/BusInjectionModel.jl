"""
    build_bim_polar!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{DC})

Build the classic DC OPF model. 

Variables:
- `pj` net bus power injection variable
- `p_gen` all generators are given real power injection decision variables at their busses.
- `v_ang` voltage angle at all busses
"""
function build_bim_polar!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{DC})
    T = net.Ntimesteps
    Y, trmnls = Ysparse(net)
    # this is a single phase model so we only need to ensure that we index variables by bus in the
    # same order as the admittance matrix
    y_busses = [trm.bus for trm in trmnls]
    B = imag(Y)
    gen_busses = generator_busses(net)
    batt_busses = battery_busses(net)

    @variable(m, pj[y_busses, 1:T])
    @variable(m, v_ang[y_busses, 1:T])

    # busses can have 1 to many Generators and or Batteries
    # so we make variables with indices (Bus, String, Time)
    m[:p_gen] = Dict{String, Dict{String, Vector{JuMP.VariableRef}}}()
    if !isempty(gen_busses)
        for b in gen_busses
            m[:p_gen][b] = Dict{String, Vector{JuMP.VariableRef}}()
            for (i, gen) in enumerate(net[b][:Generator])
                key = ismissing(gen.name) ? string(i) : gen.name
                m[:p_gen][b][key] = @variable(m, [1:T], lower_bound = gen.pmin, upper_bound = gen.pmax)
            end
        end
        # document the generator bus power injection variable
        # TODO a CommonOPF.Dimension for resources? The last dimension here is an integer for the
        # generator index but we could use a Generator.name if available
        net.var_info[:p_gen] = CommonOPF.VariableInfo(
            :p_gen,
            "generator bus power injection ",
            CommonOPF.RealPowerUnit,
            (CommonOPF.BusDimension, CommonOPF.TimeDimension)
        )
    end
    m[:p_discharge] = Dict{String, Dict{String, Vector{JuMP.VariableRef}}}()
    m[:p_charge] = Dict{String, Dict{String, Vector{JuMP.VariableRef}}}()
    m[:soc] = Dict{String, Dict{String, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 1}}}()
    if !isempty(batt_busses)
        for b in batt_busses
            m[:p_discharge][b] = Dict{String, Vector{JuMP.VariableRef}}()
            m[:p_charge][b] = Dict{String, Vector{JuMP.VariableRef}}()
            m[:soc][b] = Dict{String, JuMP.Containers.DenseAxisArray{JuMP.VariableRef, 1}}()
            # TODO rather than for loop define constraints using a set of batteries if possible
            for (i, batt) in enumerate(net[b][:Battery])
                key = ismissing(batt.name) ? string(i) : batt.name
                m[:p_discharge][b][key] = @variable(
                    m, [1:T], lower_bound = 0, upper_bound = batt.pmax / net.Sbase
                )
                m[:p_charge][b][key] = @variable(
                    m, [1:T], lower_bound = 0, upper_bound = abs(batt.pmin) / net.Sbase
                )
                m[:soc][b][key] = @variable(
                    m, [0:T], 
                    lower_bound = batt.soc_min / net.Sbase,
                    upper_bound = batt.soc_max / net.Sbase
                )

                @constraint(m, m[:soc][b][key][0] == batt.soc_init / net.Sbase)

                @constraint(m, [t = 1:T],
                    m[:soc][b][key][t] == m[:soc][b][key][t-1] + 
                        (1 / net.time_steps_per_hour) * (
                            batt.charge_efficiency * m[:p_charge][b][key][t] -
                            m[:p_discharge][b][key][t] / batt.discharge_efficiency
                        )
                )
                
                @constraint(m, [t = 1:T],
                    m[:soc][b][key][t-1] + (1 / net.time_steps_per_hour) * batt.charge_efficiency * m[:p_charge][b][key][t]
                    <= batt.soc_max / net.Sbase
                )

                @constraint(m, [t = 1:T],
                    m[:soc][b][key][t-1] - (1 / net.time_steps_per_hour) / batt.discharge_efficiency * m[:p_discharge][b][key][t]
                    >= batt.soc_min / net.Sbase
                )
            end
        end
        # document the battery charge and discharge variables
        net.var_info[:p_discharge] = CommonOPF.VariableInfo(
            :p_discharge,
            "non-negative battery discharge power injection, upper-bounded by the battery's maximum discharge rate",
            CommonOPF.RealPowerUnit,
            (CommonOPF.BusDimension, CommonOPF.TimeDimension)
        )
        net.var_info[:p_charge] = CommonOPF.VariableInfo(
            :p_charge,
            "non-negative battery charge power injection, upper-bounded by the absolute value of the battery's minimum charge rate",
            CommonOPF.RealPowerUnit,
            (CommonOPF.BusDimension, CommonOPF.TimeDimension)
        )
        @warn "Battery SOC evolution is defined via constraints but no binaries have been added to prevent simultaneous charging and discharging"
    end

    # document the variables
    net.var_info[:pj] = CommonOPF.VariableInfo(
        :pj,
        "real net bus power injection ",
        CommonOPF.RealPowerUnit,
        (CommonOPF.BusDimension, CommonOPF.TimeDimension)
    )
    net.var_info[:v_ang] = CommonOPF.VariableInfo(
        :v_ang,
        "voltage angle",
        CommonOPF.RadiansUnit,
        (CommonOPF.BusDimension, CommonOPF.TimeDimension)
    )

    net.constraint_info[:reference_angle] = CommonOPF.ConstraintInfo(
        :reference_angle,
        "voltage angle at the substation bus fixed to zero",
        MOI.EqualTo{Float64},
        (CommonOPF.TimeDimension,)
    )
    net.constraint_info[:load_balance] = CommonOPF.ConstraintInfo(
        :load_balance,
        "total net real power injection in the network equals zero",
        MOI.EqualTo{Float64},
        (CommonOPF.TimeDimension,)
    )
    net.constraint_info[:bus_real_power_injection_constraints] = CommonOPF.ConstraintInfo(
        :bus_real_power_injection_constraints,
        "net real power injection at each bus equals the negative susceptance matrix times voltage angles",
        MOI.EqualTo{Float64},
        (CommonOPF.TimeDimension, CommonOPF.BusDimension)
    )

    # reference angle set to zero
    m[:reference_angle] = @constraint(m, v_ang[net.substation_bus, :] .== 0.0)

    # net power in network must be zero (load balance)
    m[:load_balance] = @constraint(m, [t in 1:T], sum(pj[:, t]) == 0.0)

    m[:bus_real_power_injection_constraints] = @constraint(
        m, [t in 1:T],
        pj[:, t] .== -B * v_ang[:,t]
    )

    # set the power injections to loads (can be zero) and include generators as applicable
    m[:pj_constraints] = Dict()
    for j in y_busses

        sj = sj_per_unit(j, net)  # load
        # TODO slack for loads
        p_load = real(sj)

        p_gen = zeros(T)
        if j in gen_busses
            # generator dispatch decisions summed across the bus
            p_gen = sum(values(m[:p_gen][j]))
        end
        p_discharge = zeros(T)
        p_charge = zeros(T)
        if j in batt_busses
            # battery discharge decisions summed across the bus
            p_discharge = sum(values(m[:p_discharge][j]))
            # battery charge decisions summed across the bus
            p_charge = sum(values(m[:p_charge][j]))
        end

        # order of this constraint determines sign of JuMP.shadow_price, beware
        m[:pj_constraints][j] = @constraint(m, [t in 1:T],
            0.0 == pj[j,t] - p_load[t] - p_gen[t] - p_discharge[t] + p_charge[t]
        )

    end

    m[:line_limits] = Dict()
    m[:line_flow] = Dict()
    # fill in line_flow expressions and limit them as applicable
    for e in edges(net)
        i,j = indexin(e, y_busses)

        if net[e] isa CommonOPF.Conductor

            m[:line_flow][e] = @expression(m, [t in 1:T], B[i,j] * (v_ang[e[1], t] - v_ang[e[2], t]))

            if !ismissing(net[e].amps_limit)
                m[:line_limits][e] = @constraint(m, [t in 1:T],
                    -net[e].amps_limit <=  m[:line_flow][e][t] <= net[e].amps_limit
                )
            end

        elseif net[e] isa CommonOPF.ParallelConductor
            # we apply individual conductor limits as available
            # NOTE that the total susceptance of the conductors in edge `e` is accounted for in
            # matrix `B` above
            m[:line_flow][e] = Dict()
            m[:line_limits][e] = []

            for (idx, cond) in enumerate(net[e].conductors)
                bij = CommonOPF.susceptance(cond, CommonOPF.SinglePhase)

                # we use cond.name as a key in m[:line_flow][e] if the name is not missing
                if !(ismissing(cond.name))
                    idx = cond.name
                end

                m[:line_flow][e][idx] = @expression(m, [t in 1:T],
                    bij * (v_ang[e[1], t] - v_ang[e[2], t])
                )

                if !ismissing(cond.amps_limit)
                    push!(m[:line_limits][e],
                        @constraint(m, [t in 1:T],
                            -cond.amps_limit <= m[:line_flow][e][idx][t] <= cond.amps_limit
                        )
                    )
                end
            end
        end
    end
    nothing
end
