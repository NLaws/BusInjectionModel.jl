

# TODO injection_at_bus is in BranchFlowModel.jl too; move it to CommonOPF and put _per_unit in name
"""
    injection_at_bus(j::String, net::Network{MultiPhase})::Tuple{Vector{Vector{<:Real}}, Vector{Vector{<:Real}}}

return the real and reactive power injections as vectors with 3 phase indices and net.Ntimesteps time
indices like:
```julia
Pj, Qj = injection_at_bus(my_bus, net)
...
Pj[phase][time_step]
```
"""
function injection_at_bus(j::String, net::Network{MultiPhase})::Tuple{Vector{Vector{<:Real}}, Vector{Vector{<:Real}}}
    Pj = [zeros(net.Ntimesteps) for _ in 1:3] # first dim is phase, like Pj[phs][t]
    Qj = [zeros(net.Ntimesteps) for _ in 1:3]
    if j in real_load_busses(net)
        for phs in 1:3
            Pj[phs] = -net[j, :kws, phs] * 1e3 / net.Sbase
        end
    end
    if j in reactive_load_busses(net)
        for phs in 1:3 
            Qj[phs] = -net[j, :kvars, phs] * 1e3 / net.Sbase
        end
    end
    return Pj, Qj
end


# TODO substation_voltage is in BranchFlowModel.jl too; move it to CommonOPF 
"""
    substation_voltage(net::Network{MultiPhase})::Vector{ComplexF64}

Parse `net.v0` into a Vector{ComplexF64}, allowing for `net.v0` to be a `Real`,
`AbstractVector{<:Real}`, or `AbstractVector{<:Complex}`.
"""
function substation_voltage(net::Network{MultiPhase})::Vector{ComplexF64}

    if typeof(net.v0) <: Real
        return [
            net.v0 + 0im; 
            -0.5*net.v0 - im*sqrt(3)/2 * net.v0; 
            -0.5*net.v0 + im*sqrt(3)/2 * net.v0
        ]

    elseif typeof(net.v0) <: AbstractVector{<:Real}
        return [
            net.v0[1] + 0im; 
            -0.5 * net.v0[2] - im*sqrt(3)/2 * net.v0[2]; 
            -0.5 * net.v0[3] + im*sqrt(3)/2 * net.v0[3]
        ]

    elseif typeof(net.v0) <: AbstractVector{<:Complex}
        return net.v0

    else  
        throw(@error "unsupported type for Network.v0 $(typeof(net.v0))")
    end

    return w0
end


# TODO phi_ij is in BranchFlowModel.jl too; move it to CommonOPF 
"""
    phi_ij(j::String, net::Network, v::AbstractVector)

Down-select the vector v by the phase from i -> j
"""
function phi_ij(j::String, net::Network, v::AbstractVector)
    n = convert(Vector{GenericAffExpr{ComplexF64, VariableRef}}, [0im; 0im; 0im])
    for x in phases_into_bus(net, j)
        n[x] = v[x]
    end
    return n
end


# TODO cj is in BranchFlowModel.jl too; move it to CommonOPF 
"""
    cj(A)

short cut for conj(transpose(A))
"""
function cj(A)
    conj(transpose(A))
end


# TODO phases_of_vector  move it to CommonOPF 
"""
    phases_of_vector(v::AbstractVector{T}, phases::AbstractVector{Int}) where T

Return only the `phases` of vector `v`.
"""
function phases_of_vector(v::AbstractVector{T}, phases::AbstractVector{Int}) where T
    v2 = T[]
    for i in phases
        push!(v2,v[i])
    end
    return v2
end


"""
    build_bim_rectangular!(m::JuMP.AbstractModel, net::Network{MultiPhase}, ::Val{Unrelaxed})

Model builder for multi-phase, unrelaxed BIM with rectangular voltage variables. See the 
    [Multi-Phase Bus Injection Model (Unrelaxed)](@ref) math for details.

Adds the variables:
- `m[:v]` with complex values for all busses in `CommonOPF.busses(net)`
- `m[:s0]` for the complex slack bus power injection
"""
function build_bim_rectangular!(m::JuMP.AbstractModel, net::Network{MultiPhase}, ::Val{Unrelaxed})
    Tr = phases_of_vector  # "Tr" for transform
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

        Pj, Qj = injection_at_bus(j, net)
        Sj = Pj + im * Qj
        Sj = hcat(Sj...)  # time X phase

        phases = phases_connected_to_bus(net, j)
        # have to down-select to phases to avoid bad constraints like 0 = (1.7 - 0.8im)
        @constraint(m, [t in 1:T],
            Tr(Sj[t, :], phases) .== Tr(sum(
                diag(
                    v[j][t] * cj( v[j][t]  - phi_ij(j, net, v[k][t]) ) * cj(yij_per_unit(j, k, net))
                )
                for k in connected_busses(j, net)
            ), phases)
        )    
    end

end
