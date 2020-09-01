# Description

Yarf is an open-source circuit simulator written in Python. The package contains an API for netlist description and to create analyses to run in the circuit. One goal of the project is to integrate Yarf with scikit-rf and openEMS to create a powerful open-source development environment for RF engineers.

# Dependencies

`conda install python numpy scipy matplotlib`

# TODO List

## Analysis Related:
[ ] transient analysis (trapezoidal, variable time-step, etc)
[ ] fft for transient results
[ ] harmonic balance single-tone and multitone (and for autonomous circuits)
[ ] PSS analysis (shooting)
[ ] PZ analysis
[ ] S-parameter and LSSP
[ ] AC noise analysis
[ ] transient noise sources
[ ] nodeset for DC analysis
[x] AC analysis
[x] gmin and source stepping continuation methods
[x] non-linear DC analysis
[x] linear DC analysis
[x] creation of MNA matrices using element stamps

## Modelling Related:
[ ] add ac model of controlled sources
[ ] fix/differentiate dc and ac model for vsource and isources
[ ] implement mosfet (level 3 model and maybe some bsim boy, see hspice doc)
[ ] implement transformer
[ ] implement current probe
[ ] transmission lines
[ ] S-Parameter model (AC analysis and ChirpZ for transient)
[ ] noise and temperature modelling
[x] AC models for AC analysis
[x] implement diode model
[x] implement bjt spice gummel-poon (TODO: check for singularities)

## API and Code Related:
[ ] need to add exception handling instead of all the if and elses
[ ] send solve_linear and solve_dc_nonlinear to a Solver.py file so it can be reused
[ ] add check at the top of each Yarf.add_XX() to see if the device or analysis name already exists in the netlist
[ ] add a method to check if the netlist has at least one purposely placed gnd node
[ ] implement a test suite comparing the results with Xyce
[ ] add subcircuit support
[ ] add a remove method for devices and analyses
[ ] make the netlist not case sensitive (? not sure)
[x] declare logger objects someplace else
[x] add easy methods to access the results (e.g. get_v('R1'))
[x] create git repository
[x] API for netlist creation
[x] move element stamps from DC.py to device classes
[x] add logging support (for debugging)


