# The fixed point magnitude and phasor models
# v = −Y_LL^−1 Y_L0 v_0 + Y_LL^−1 diag(v_FP^*) s^*

# TODO quadratic model with amperage variables and limits (voltage rectangular phasor)

# Can index JuMP variable on BusTerminal! but should I?


"""
    add_complex_terminal_voltage_variable(
        m::JuMP.Model, 
        net::Network{MultiPhase}, 
        trmnls::Vector{CommonOPF.BusTerminal},
    )

Add a variable with symbol `:v` for the complex voltage indexed on `trmnls` and `1:net.Ntimesteps`
"""
function add_complex_terminal_voltage_variable(
        m::JuMP.Model, 
        net::Network{MultiPhase}, 
        trmnls::Vector{CommonOPF.BusTerminal},
    )
    @variable(m, v[trmnls, 1:net.Ntimesteps], set=ComplexPlane())

    net.var_info[:v] = CommonOPF.VariableInfo(
        :v,
        "complex terminal voltage",
        CommonOPF.VoltUnit,
        (CommonOPF.BusTerminalDimension, CommonOPF.TimeDimension)
    )
    return v
end


function add_complex_terminal_power_variable(m::JuMP.Model, net::Network{MultiPhase})
    ts = terminals(net)
    @variable(m, s[ts, 1:net.Ntimesteps], set=ComplexPlane())

    net.var_info[:s] = CommonOPF.VariableInfo(
        :s,
        "complex terminal net power injection",
        CommonOPF.ComplexPowerUnit,
        (CommonOPF.BusTerminalDimension, CommonOPF.TimeDimension)
    )
end


"""
     add_or_update_fixed_point_constraint(
        m::JuMP.Model, 
        net::Network{MultiPhase}, 
        v_fp::Matrix{ComplexF64},
    )

Apply the fixed point voltage constraint to the model `m`. If `:fixed_point_con in keys(m.obj_dict)`
then the existing constraints are deleted first.
"""
function add_or_update_fixed_point_constraint(
        m::JuMP.Model, 
        net::Network{MultiPhase}, 
        v_fp::Matrix{ComplexF64},
    )

    v0 = substation_voltage(net)
    
    # we use the known loads as the fixed point power
    s_fp = terminals_sj_per_unit(net, m[:y_terminals])[m[:ll_indices]]
    s_fp = hcat(s_fp...)

    # TODO s variable
    if :fixed_point_con in keys(m.obj_dict)
        for term in m[:ll_terminals], t in 1:net.Ntimesteps
            # TODO? should not have matrix in last dimension?
            JuMP.delete(m, m[:fixed_point_con][t][term, 1])
        end
        delete!(m.obj_dict, :fixed_point_con)
    end
    @constraint(m, fixed_point_con[t in 1:net.Ntimesteps],
        m[:v][:, t] .== -m[:inv_Yll] * m[:Y_l0] * v0 .+ m[:inv_Yll] * diagm(conj(v_fp[:, t]))^-1 * conj(s_fp[t, :])
    )
    # document the constraints
    # TODO update dimensions once the todo about the 1 in the last dimension is addressed
    c = m[:fixed_point_con][1][m[:ll_terminals][1], 1]  # time step 1
    net.constraint_info[:bus_power_injection_constraints] = CommonOPF.ConstraintInfo(
        :fixed_point_con,
        "fixed point voltage equation",
        typeof(MOI.get(m, MOI.ConstraintSet(), c)),
        (CommonOPF.TimeDimension, CommonOPF.BusTerminalDimension),
    )

    nothing
end


"""
    build_bim_rectangular!(m::JuMP.Model, net::Network{MultiPhase}, ::Val{FixedPointLinear})

Build the fixed point linear model.
"""
function build_bim_rectangular!(m::JuMP.Model, net::Network{MultiPhase}, ::Val{FixedPointLinear})
    Y, y_terminals = Ysparse(net)
    # N = size(Y, 1)
    m[:y_terminals] = y_terminals
    source_bus_terminals = filter(x -> x.bus == net.substation_bus, y_terminals)
    source_indices = [term.Y_index for term in source_bus_terminals]

    ll_terminals = m[:ll_terminals] = filter(x -> x.bus != net.substation_bus, y_terminals)
    ll_indices = m[:ll_indices] = [term.Y_index for term in ll_terminals]

    Y_ll = Y[ll_indices, ll_indices]
    m[:Y_l0] = Y[ll_indices, source_indices] * net.Zbase
    v0 = substation_voltage(net)

    # TODO way to not make dense matrix?
    m[:inv_Yll] = inv(Matrix(Y_ll)) / net.Zbase

    N_ll = length(ll_indices)
    v_fp = ones(N_ll, net.Ntimesteps) * 0im
    for (i, term) in enumerate(ll_terminals), t in 1:net.Ntimesteps
        v_fp[i, t] = v0[term.phase]
    end
    
    # we use the known loads as the fixed point power
    s_fp = terminals_sj_per_unit(net, y_terminals)[ll_indices]
    s_fp = hcat(s_fp...)

    v = add_complex_terminal_voltage_variable(m, net, ll_terminals)
    # initialize voltage with v0
    for term in ll_terminals, t in 1:net.Ntimesteps
        set_start_value(real(v[term, t]), real(v0[term.phase]))
        set_start_value(imag(v[term, t]), imag(v0[term.phase]))
    end
    
    # TODO s variable
    add_or_update_fixed_point_constraint(m, net, v_fp)

    nothing
end


function max_abs_diff(v1::AbstractMatrix, v2::AbstractMatrix)
    maximum(abs.(v1 .- v2))
end


function solve_fixed_point_to_tol(m::JuMP.Model, net::Network{MultiPhase}, tol::Float64=1e-3, max_iter::Int=10)

    vdiff = tol + 1.0
    iter = 1

    # get the first voltage
    optimize!(m)
    v_prev = value.(m[:v])

    # perform at least one iteration to get vdiff
    while iter <= max_iter && vdiff > tol
        add_or_update_fixed_point_constraint(m, net, v_prev.data)
        optimize!(m)
        v_now = value.(m[:v])

        vdiff = max_abs_diff(v_prev.data, v_now.data)
        iter += 1
        v_prev = v_now
        @debug("iter $iter \t max_abs_diff $(@sprintf("%.6f", vdiff))")
    end

    if vdiff > tol
        @warn "Model not solved to tolerance of $tol. vdiff is $vdiff"

    else
        @info "Model solved to tolerance of $tol in $(iter-1) iterations."
    end

    nothing
end
