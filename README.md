# Description

YalRF is an open-source circuit simulator written in Python. The package contains an API for netlist description and circuit simulation. The generated data can be easily post-processed using numpy and Jupyter Notebooks.

The main goal of this project is to implement a stable and powerful multi-tone harmonic balance engine with support to autonomous circuits.

Another goal will be to integrate YalRF with the scikit-rf / openEMS / SignalIntegrity packages to create a powerful open-source development environment for RF engineers.

Example of usage:
```python
from yalrf import YalRF

y = YalRF('Voltage Divider')

y.add_resistor('R1', 'n1', 'n2', 100)
y.add_resistor('R2', 'n2', 'gnd', 25)

v1 = y.add_vdc('V1', 'n1', 'gnd', 1)

dc = y.add_dc_analysis('DC1')

y.run('DC1')
y.print_dc_voltages('DC1')
y.print_dc_currents('DC1')
```

# Dependencies

`conda install python numpy scipy matplotlib`

# TODO List

## Analysis Related:
- [ ] transient analysis (trapezoidal, variable time-step, initial condition (??), etc)
- [ ] fft for transient results
- [ ] harmonic balance single-tone and multitone (and for autonomous circuits)
- [ ] PSS analysis (shooting)
- [ ] PZ analysis
- [ ] S-parameter and LSSP
- [ ] AC noise analysis
- [ ] transient noise sources
- [ ] nodeset for DC analysis
- [x] AC analysis
- [x] gmin and source stepping continuation methods
- [x] non-linear DC analysis
- [x] linear DC analysis
- [x] creation of MNA matrices using element stamps

## Modelling Related:
- [ ] implement mosfet
- [ ] implement transformer
- [ ] implement current probe
- [ ] implement ideal operational amplifier
- [ ] implement controlled switch (relay)
- [ ] transmission lines
- [ ] S-Parameter block (AC analysis and ChirpZ for transient)
- [ ] noise and temperature modelling
- [x] add ac model of controlled sources
- [x] AC models for AC analysis
- [x] implement diode model
- [x] implement bjt spice gummel-poon (TODO: check for singularities)

## API and Code Related:
- [ ] current of bjts and diodes are currently stored inside the model (need to move to the analysis)
- [ ] review log messages to see if the type make sense (info, warning, error)
- [ ] exception handling instead of all the ifs and elses
- [ ] send solve_linear and solve_dc_nonlinear to a Solver.py file so it can be reused
- [ ] check at the top of each YalRF.add_device() call if the device name already exists in the netlist
- [ ] check if the netlist has at least one purposely placed gnd node
- [ ] subcircuit support
- [ ] remove method for devices and analyses in YalRF API
- [ ] input sanitizaiton can be largely improved
- [ ] there is a lot of waste in memory and performance that can be improved
- [x] implement a test suite comparing the results with Xyce
- [x] fix/differentiate dc and ac model for vsource and isource
- [x] declare logger objects someplace else
- [x] add easy methods to access the results (e.g. get_v('R1'))
- [x] create git repository
- [x] API for netlist creation
- [x] move element stamps from DC.py to device classes
- [x] add logging support (for debugging)


