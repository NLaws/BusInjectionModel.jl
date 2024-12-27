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

    # s_j = sum_{k: j~k} Y[j,k]^* ( |v_j|^2 - v_j v_k^*)
    for j in busses(net)

        if j == net.substation_bus

            @constraint(m, [t in 1:T],
                s0[j][t] .== sum(
                    diag(
                        v[j][t] * cj( v[j][t]  - phi_ij(j, net, v[k][t]) ) * cj(yij_per_unit(j, k, net))
                    )
                    for k in j_to_k(j, net)
                )
            )
            continue
        end

        Pj, Qj = sj_per_unit(j, net)
        Sj = Pj + im * Qj
        Sj = hcat(Sj...)  # time X phase

        phases = phases_connected_to_bus(net, j)
        # have to down-select to phases to avoid bad constraints like 0 = (1.7 - 0.8im)
        @constraint(m, [t in 1:T],
            Sj[t, phases] .== sum(
                diag(
                    v[j][t] * cj( v[j][t]  - phi_ij(j, net, v[k][t]) ) * cj(yij_per_unit(j, k, net))
                )[phases]
                for k in connected_busses(j, net)
            )
        )    
    end

end
