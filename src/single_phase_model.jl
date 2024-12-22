"""
    build_bim!(m::JuMP.AbstractModel, net::Network, mtype::ModelType=Unrelaxed)

Top-level model builder that dispatches the ModelType enum
"""
function build_bim!(m::JuMP.AbstractModel, net::Network, mtype::ModelType=Unrelaxed)
    build_bim!(m::JuMP.AbstractModel, net::Network, Val(mtype))
end


"""
    build_bim!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{Unrelaxed})

Model builder for single-phase, unrelaxed BIM with rectangular voltage variables. See the 
    [Single Phase Bus Injection Model (Unrelaxed)](@ref) math for details.
"""
function build_bim!(m::JuMP.AbstractModel, net::Network{SinglePhase}, ::Val{Unrelaxed})
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
