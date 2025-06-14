# Single Phase Models
BusInjectionModel.jl provides methods to build many different variations of the Branch Flow Model,
including single phase and multiphase models. Each of the single phase model types supported are documented below.
```@contents
Pages = ["single_phase_models.md"]
Depth = 3
```
```@setup imports
using BusInjectionModel
using CommonOPF
using JuMP
```


## `Unrelaxed` models
The `Unrelaxed` multiphase model is built by passing a `JuMP.Model`, `Network{SinglePhase}`, and the
`Unrelaxed` type to [`build_bim_rectangular!`](@ref) or [`build_bim_polar!`](@ref).


### Rectangular Voltage
```@example imports
net = CommonOPF.Network_IEEE13_SinglePhase()
m = JuMP.Model()

build_bim_rectangular!(m, net, Unrelaxed)
println("Variable information:")
CommonOPF.print_var_info(net)
println("Constraint information:")
CommonOPF.print_constraint_info(net)
```

The math underlying the model is as follows:
```math
s_j = \sum_{k: j \sim k} Y_{jk}^* \left( |v_j|^2 - v_j v_k^* \right)
\quad \forall j \in \mathcal{N}
```
For the nomenclature see TODO.


### Polar Voltage
```@example imports
net = CommonOPF.Network_IEEE13_SinglePhase()
m = JuMP.Model()

build_bim_polar!(m, net, Unrelaxed)
println("Variable information:")
CommonOPF.print_var_info(net)
println("Constraint information:")
CommonOPF.print_constraint_info(net)
```

The math underlying the model is as follows:
```math
\begin{aligned}
&p_j =  |v_j| \sum_{i \in 1\dots\mathcal{N}} |v_i| \left[
    G_{\left[i,j\right]} \cos(\angle v_j - \angle v_i) + B_{\left[i,j\right]} \sin(\angle v_j - \angle v_i)
    \right]
    \quad \forall j \in \mathcal{N}
\\
&q_j =  |v_j| \sum_{i \in 1\dots\mathcal{N}} |v_i| \left[
    G_{\left[i,j\right]} \sin(\angle v_j - \angle v_i) - B_{\left[i,j\right]} \cos(\angle v_j - \angle v_i)
    \right]
    \quad \forall j \in \mathcal{N}
\\
&|v_{\text{substation bus}}| = v_0 
\\
&\angle v_{\text{substation bus}} = 0
\end{aligned}
```