"""
    build_bim!(m::JuMP.AbstractModel, net::Network{MultiPhase}, mtype::ModelType=FixedPointLinear)

Top-level multiphase builder that dispatches the ModelType enum
"""
function build_bim!(m::JuMP.AbstractModel, net::Network{MultiPhase}, mtype::ModelType=FixedPointLinear)
    build_bim!(m::JuMP.AbstractModel, net::Network{MultiPhase}, Val(mtype))
end
