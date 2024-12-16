# BusInjectionModel.jl

Documentation for BusInjectionModel.jl

A work-in-progress package for optimal power flow (OPF) models based on the bus injection model.

This package is part of a federation of packages to support OPF modeling:
- [CommonOPF](https://github.com/NLaws/CommonOPF.jl) provides the basis for BusInjectionModel.jl
  and:
- [BranchFlowModel](https://github.com/NLaws/BranchFlowModel.jl), which provides a similar interface
  to BusInjectionModel.jl.

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
