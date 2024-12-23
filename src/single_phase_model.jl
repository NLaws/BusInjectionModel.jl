"""
    build_bim_rectangular!(m::JuMP.AbstractModel, net::Network, mtype::ModelType=Unrelaxed)

Top-level model builder that dispatches the ModelType enum
"""
function build_bim_rectangular!(m::JuMP.AbstractModel, net::Network, mtype::ModelType=Unrelaxed)
    build_bim_rectangular!(m::JuMP.AbstractModel, net::Network, Val(mtype))
end


"""
    build_bim_rectangular!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{Unrelaxed})

Model builder for single-phase, unrelaxed BIM with rectangular voltage variables. See the 
    [Single Phase Bus Injection Model (Unrelaxed)](@ref) math for details.

Adds the variables:
- `m[:v]` with complex values for all busses in `CommonOPF.busses(net)`
- `m[:s0]` for the complex slack bus power injection
"""
function build_bim_rectangular!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{Unrelaxed})
    T = net.Ntimesteps

    add_time_vector_variables!(m, net, :v, busses(net); set=ComplexPlane)
    v = m[:v]
    @variable(m, set = ComplexPlane(), s0[1:T])

    # s_j = sum_{k: j~k} Y[j,k]^* ( |v_j|^2 - v_j v_k^*)
    for j in busses(net)

        if j == net.substation_bus

            @constraint(m, [t in 1:T],
                s0[t] == sum(
                    conj(yij_per_unit(j, k, net)) * (
                        real( v[j][t] )^2 + imag( v[j][t] )^2 - v[j][t] * conj(v[k][t])
                    )
                    for k in j_to_k(j, net)
                )
            )

            @constraint(m, [t in 1:T],
                v[j][t] == net.v0[t]
            )
            continue
        end

        @constraint(m, [t in 1:T],
            sj_per_unit(j, net)[t] == sum(
                conj(yij_per_unit(j, k, net)) * (
                    real( v[j][t] )^2 + imag( v[j][t] )^2 - v[j][t] * conj(v[k][t])
                )
                for k in connected_busses(j, net)
            )
        )
    
    end

end


"""
    build_bim_polar!(m::JuMP.AbstractModel, net::Network, mtype::ModelType=Unrelaxed)

Top-level model builder that dispatches the ModelType enum
"""
function build_bim_polar!(m::JuMP.AbstractModel, net::Network, mtype::ModelType=Unrelaxed)
    build_bim_polar!(m::JuMP.AbstractModel, net::Network, Val(mtype))
end


"""
    build_bim_polar!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{Unrelaxed})


"""
function build_bim_polar!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{Unrelaxed})
    T = net.Ntimesteps
    add_time_vector_variables!(m, net, :v_mag, busses(net))
    add_time_vector_variables!(m, net, :v_ang, busses(net))
    v_mag = m[:v_mag]
    v_ang = m[:v_ang]

    # voltage angles start at zero and between -π and π
    for b in setdiff(busses(net), [net.substation_bus])
        for t in 1:T
            JuMP.set_start_value(v_ang[b][t], 0.0)
            JuMP.set_lower_bound(v_ang[b][t], -π)
            JuMP.set_upper_bound(v_ang[b][t], π)
            if !(b in generator_busses(net))
                JuMP.set_start_value(v_mag[b][t], 1.0)
                JuMP.set_lower_bound(v_mag[b][t], 0.9)
                JuMP.set_upper_bound(v_mag[b][t], 1.1)
            end
        end
    end

    # slack bus voltage and variables
    @constraint(m, [t in 1:T], v_mag[net.substation_bus][t] == net.v0)
    @constraint(m, [t in 1:T], v_ang[net.substation_bus][t] == 0.0)
    @variable(m, p0[1:T])
    @variable(m, q0[1:T])

    # busses with known power injections
    @constraint(m, con_real_power[j in real_load_busses(net), t in 1:T],
        real(sj_per_unit(j, net)[t]) == (
            v_mag[j][t]
            * sum( v_mag[i][t] * (
                    real(Yij_per_unit(i, j, net)) # conductance
                *   cos(v_ang[j][t] - v_ang[i][t])
                +   imag(Yij_per_unit(i, j, net)) # susceptance
                *   sin(v_ang[j][t] - v_ang[i][t])
            ) for i in union(connected_busses(j, net), [j])
            )
        )
    )

    @constraint(m, con_reactive_power[j in reactive_load_busses(net), t in 1:T],
        imag(sj_per_unit(j, net)[t]) == (
            v_mag[j][t] 
            * sum( v_mag[i][t] * (
                    real(Yij_per_unit(i, j, net)) # conductance
                *   sin(v_ang[j][t] - v_ang[i][t])
                -   imag(Yij_per_unit(i, j, net)) # susceptance
                *   cos(v_ang[j][t] - v_ang[i][t])
            ) for i in union(connected_busses(j, net), [j])
            )
        )
    )

    #P-V bus Q variables and voltage constraints
    if !isempty(generator_busses(net))
        add_time_vector_variables!(m, net, :q_gen, generator_busses(net))

        for j in generator_busses(net)

            @constraint(m, [t in 1:T],
                net[j][:Generator].kws1[t] * 1e3 / net.Sbase == (
                    v_mag[j][t]
                    * sum( v_mag[i][t] * (
                            real(Yij_per_unit(i, j, net)) # conductance
                        *   cos(v_ang[j][t] - v_ang[i][t])
                        +   imag(Yij_per_unit(i, j, net)) # susceptance
                        *   sin(v_ang[j][t] - v_ang[i][t])
                    ) for i in union(connected_busses(j, net), [j])
                    )
                )
            )

            @constraint(m, [t in 1:T],
                m[:q_gen][j][t] == (
                    v_mag[j][t]
                    * sum( v_mag[i][t] * (
                            real(Yij_per_unit(i, j, net)) # conductance
                        *   sin(v_ang[j][t] - v_ang[i][t])
                        -   imag(Yij_per_unit(i, j, net)) # susceptance
                        *   cos(v_ang[j][t] - v_ang[i][t])
                    ) for i in union(connected_busses(j, net), [j])
                    )
                )
            )

            @constraint(m, [t in 1:T], v_mag[j][t] == net[j][:Generator].voltage_pu[t])

        end
    end


    nothing
end
