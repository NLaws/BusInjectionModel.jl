# BusInjectionModel.jl

Documentation for BusInjectionModel.jl

A work-in-progress package for optimal power flow (OPF) models based on the bus injection model.

This package is part of a federation of packages to support OPF modeling:
- [CommonOPF](https://github.com/NLaws/CommonOPF.jl) provides the basis for BusInjectionModel.jl
  and:
- [BranchFlowModel](https://github.com/NLaws/BranchFlowModel.jl), which provides a similar interface
  to BusInjectionModel.jl.

# Inputs
Inputs are defined using [`CommonOPF.Network` structs](https://nlaws.github.io/CommonOPF.jl/stable/network/). 


# Building a Model
Building a `BranchFlBusInjectionModelowModel` requires three things:
1. a JuMP Model,
2. a `CommonOPF.Network`, and
3. the type of model to be built, i.e. one of the [`BusInjectionModel.ModelType`](@ref)
```@docs
BusInjectionModel.ModelType
```
To build a model see [`build_bim_rectangular!`](@ref) and [`build_bim_polar!`](@ref)

