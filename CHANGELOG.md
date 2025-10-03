# BusInjectionModel Change log

## v0.1.2
- add fixed point linear model, multi and single phase
    - leverage `CommonOPF.BusTerminal` as model indices s.t. do not have to implement single and
      multiphase model builders (requires CommonOPF v0.5.0)

## v0.1.1
- add multiphase rectangular model

## v0.1.0 First release 
- single phase models for rectangular and polar voltages