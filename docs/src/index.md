# BusInjectionModel.jl

Documentation for BusInjectionModel.jl


```julia
using BusInjectionModel
using CommonOPF
using JuMP
using Ipopt

m = JuMP.Model(Ipopt.Optimizer)
net = CommonOPF.Network("path/to/network/yaml-or-opendss-file")
BusInjectionModel.build_bim!(m, net, Unrelaxed)
# set objective
optimize!(m)
```
