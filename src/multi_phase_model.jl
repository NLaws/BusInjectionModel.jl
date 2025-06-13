"""
    build_bim_rectangular!(m::JuMP.AbstractModel, net::Network{MultiPhase}, ::Val{Unrelaxed})

Model builder for multi-phase, unrelaxed BIM with rectangular voltage variables. See the 
    [Multi-Phase Bus Injection Model (Unrelaxed)](@ref) math for details.

Adds the variables:
- `m[:v]` with complex values for all busses in `CommonOPF.busses(net)`
- `m[:s0]` for the complex slack bus power injection
"""
function build_bim_rectangular!(m::JuMP.AbstractModel, net::Network{MultiPhase}, ::Val{Unrelaxed})
    T = net.Ntimesteps
    v = m[:v] = multiphase_bus_variable_container()
    # bus net power injection vectors
    s0 = m[:s0] = multiphase_bus_variable_container()
    net.var_names = [:v, :s0]
    sub_v = substation_voltage(net)

    for t in 1:net.Ntimesteps
        for j in busses(net)
    
            if j == net.substation_bus
                # slack bus injection
                add_complex_vector_of_phase_variable!(
                    m, net, j, :s0, t;
                )
                # TODO s0 start value
                # slack bus voltage
                v[j][t] = sub_v
            else
                # voltage variables
                add_complex_vector_of_phase_variable!(
                    m, net, j, :v, t;
                    upper_bound_mag = net.bounds.v_upper_mag,
                    lower_bound_mag = net.bounds.v_lower_mag,
                )
                for phs in phases_connected_to_bus(net, j)
                    set_start_value(real(m[:v][j][t][phs]), real(sub_v[phs]))
                    set_start_value(imag(m[:v][j][t][phs]), imag(sub_v[phs]))
                end
            end

        end
    end

    # document the variables
    net.var_info[:v] = CommonOPF.VariableInfo(
        :v,
        "complex voltage vector",
        CommonOPF.VoltUnit,
        (CommonOPF.BusDimension, CommonOPF.TimeDimension, CommonOPF.PhaseDimension)
    )

    net.var_info[:s0] = CommonOPF.VariableInfo(
        :s0,
        "complex net bus power injection at the net.substation_bus",
        CommonOPF.ComplexPowerUnit,
        (CommonOPF.BusDimension, CommonOPF.TimeDimension, CommonOPF.PhaseDimension)
    )

    # define power
    m[:bus_power_injection_constraints] = Dict()
    for j in busses(net)

        if j == net.substation_bus

            m[:bus_power_injection_constraints][j] = @constraint(m, [t in 1:T],
                s0[j][t] .== sum(
                    diag(
                        v[j][t] * cj( v[j][t]  - phi_ij(j, net, v[k][t]) ) * cj(yij_per_unit(j, k, net))
                    )
                    for k in j_to_k(j, net)
                )
            )
            continue
        end

        sj = sj_per_unit(j, net)
        sj = hcat(sj...)  # time X phase

        phases = phases_connected_to_bus(net, j)
        # have to down-select to phases to avoid bad constraints like 0 = (1.7 - 0.8im)
        m[:bus_power_injection_constraints][j] = @constraint(m, [t in 1:T],
            sj[t, phases] .== sum(
                diag(
                    v[j][t] * cj( v[j][t]  - phi_ij(j, net, v[k][t]) ) * cj(yij_per_unit(j, k, net))
                )[phases]
                for k in connected_busses(j, net)
            )
        )    
    end

    # document the constraints
    c = m[:bus_power_injection_constraints][net.substation_bus][1][1]  # time step 1, phase 1
    net.constraint_info[:bus_power_injection_constraints] = CommonOPF.ConstraintInfo(
        :bus_power_injection_constraints,
        "net power injection definition at each bus",
        typeof(MOI.get(m, MOI.ConstraintSet(), c)),
        (CommonOPF.BusDimension, CommonOPF.TimeDimension, CommonOPF.PhaseDimension),
    )

    nothing
end
