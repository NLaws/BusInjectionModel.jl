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

    @variable(m, pj[y_busses, 1:T])
    @variable(m, v_ang[y_busses, 1:T])

    # busses can have 1 to many Generators
    # so we make variables with indices (Bus, Int, Time)
    m[:p_gen] = Dict{String, Dict{Int, Vector{JuMP.VariableRef}}}()
    if !isempty(gen_busses)
        for b in gen_busses
            m[:p_gen][b] = Dict{Int, Vector{JuMP.VariableRef}}()
            for (i, gen) in enumerate(net[b][:Generator])
                m[:p_gen][b][i] = @variable(m, [1:T], lower_bound = gen.pmin, upper_bound = gen.pmax)
            end
        end
    end

    p = m[:pj]
    v_ang = m[:v_ang]

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

    # reference angle set to zero
    @constraint(m, v_ang[net.substation_bus, :] .== 0.0)

    # net power in network must be zero (load balance)
    m[:load_balance] = @constraint(m, [t in 1:T], sum(pj[:, t]) == 0.0)

    m[:bus_real_power_injection_constraints] = @constraint(
        m, [t in 1:T],
        p[:, t] .== -B * v_ang[:,t]
    )

    # set the power injections to loads (can be zero) and include generators as applicable
    # TODO document these constraints
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

        # order of this constraint determines sign of JuMP.shadow_price, beware
        m[:pj_constraints][j] = @constraint(m, [t in 1:T],
            0.0 == p[j,t] - p_load[t] - p_gen[t]
        )

    end

    m[:line_limits] = Dict()
    m[:line_flow] = Dict()

    for e in edges(net)
        i,j = indexin(e, y_busses)

        m[:line_flow][e] = @expression(m, [t in 1:T], B[i,j] * (v_ang[e[1], t] - v_ang[e[2], t]))

        if net[e] isa CommonOPF.Conductor && !ismissing(net[e].amps_limit)
            m[:line_limits][e] = @constraint(m, [t in 1:T],
                -net[e].amps_limit <=  m[:line_flow][e][t] <= net[e].amps_limit
            )
        elseif net[e] isa CommonOPF.ParallelConductor

            if all(ismissing((c.amps_limit) for c in net[e].conductors))
                # no line limits
                continue
            end

            if !(length(net[e].conductors) == 2)
                @warn "Not setting line constraints for $(net[e]). Only 2 parallel conductors are handled."
                continue
            end

            c1 = net[e].conductors[1]
            c2 = net[e].conductors[2]

            # TODO need a different model when the conductor impedances are not equal?
            if !(CommonOPF.susceptance(c1, CommonOPF.SinglePhase) == CommonOPF.susceptance(c2, CommonOPF.SinglePhase))
                @warn "Parallel conductors without equal susceptance are not handled"
                continue
            end

            total_limit = 0.0
            if ismissing(c1.amps_limit) && !ismissing(c2.amps_limit)
                total_limit = 2 * c2.amps_limit
            elseif !ismissing(c1.amps_limit) && ismissing(c2.amps_limit)
                total_limit = 2 * c1.amps_limit
            else
                total_limit = 2 * minimum([c1.amps_limit, c2.amps_limit])
            end
            
            m[:line_limits][e] = @constraint(m, [t in 1:T],
                -total_limit <=  m[:line_flow][e][t] <= total_limit
            )
        end
    end
    nothing
end