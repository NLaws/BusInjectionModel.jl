# Multiphase Models
BusInjectionModel.jl provides methods to build many different variations of the Branch Flow Model,
including single phase and multiphase models. Each of the multihpase model types supported are documented below.
```@contents
Pages = ["multi_phase_models.md"]
Depth = 2
```
```@setup imports
using BusInjectionModel
using CommonOPF
using JuMP
```


## `Unrelaxed` models
The `Unrelaxed` multiphase model is built by passing a `JuMP.Model`, `Network{MultiPhase}`, and the
`Unrelaxed` type to [`build_bim_rectangular!`](@ref).

```@example imports
net = CommonOPF.Network_IEEE13()
m = JuMP.Model()

build_bim_rectangular!(m, net, Unrelaxed)
println("Variable information:")
CommonOPF.print_var_info(net)
println("Constraint information:")
CommonOPF.print_constraint_info(net)
```

The math underlying the model is as follows:
```math
\boldsymbol s_j^{\Phi_{j}} = \sum_{k: j \sim k} \text{diag} \left(  
    \boldsymbol v_j^{\Phi_{jk}} \left[ \boldsymbol v_j^{\Phi_{jk}} - \boldsymbol v_k^{\Phi_{jk}}  \right]^H \boldsymbol Y_{jk}^H
\right)
\quad \forall j \in \mathcal{N}
```
For the nomenclature see TODO.


## `FixedPointLinear` models
The `FixedPointLinear` multiphase model is built by passing a `JuMP.Model`, `Network{MultiPhase}`, and the
`FixedPointLinear` type to [`build_bim_rectangular!`](@ref).

```@example imports
net = CommonOPF.Network_IEEE13()
m = JuMP.Model()

build_bim_rectangular!(m, net, FixedPointLinear)
println("Variable information:")
CommonOPF.print_var_info(net)
println("Constraint information:")
CommonOPF.print_constraint_info(net)
```

The math underlying the model is as follows:
```math
\boldsymbol v = -\boldsymbol Y_{LL}^{-1} Y_{L0} \boldsymbol v_0 + Y_{LL}^{-1} \text{diag}(\boldsymbol v_{FP}^*)^{-1} \boldsymbol s^*
```
For the nomenclature see TODO.
