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
- [ ] add variable time-step for transient analysis and other integration methods
- [ ] expand two-tone harmonic balance to multi-tone
- [ ] expand harmonic balance NR for autonomous analysis instead of using scipy optimize
- [ ] PSS analysis (shooting method)
- [ ] S-parameter and LSSP
- [ ] AC noise analysis
- [ ] transient noise sources
- [ ] nodeset for DC analysis
- [ ] expand API of each analysis to process internal data

## Modelling Related:
- [ ] complete implementation of mosfet
- [ ] add transformer
- [ ] add current probe
- [ ] add controlled switch (relay)
- [ ] transmission lines
- [ ] S-Parameter and Y-Matrix blocks
- [ ] noise modelling
- [ ] more behavioral models (sum, subtract, multiply, exponential, polynomial)

## API and Code Related:
- [ ] current of bjts and diodes are currently stored inside the model (need to fix!)
- [ ] move integration to transient analysis
- [ ] review log messages to see if the type make sense (info, warning, error)
- [ ] exception handling instead of all the ifs and elses
- [ ] check at the top of each YalRF.add_device() call if the device name already exists in the netlist
- [ ] check if the netlist has at least one purposely placed gnd node
- [ ] subcircuit support
- [ ] remove devices from netlist
- [ ] input sanitizaiton can be largely improved
- [ ] there is a lot of waste in memory and performance that can be improved
- [ ] increment the test suite comparing the results with Xyce
