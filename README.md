# Description

YalRF is an open-source circuit simulator written in Python. The package contains an API for netlist description and circuit simulation. The generated data can be easily post-processed using numpy and Jupyter Notebooks.

The main goal of this project is to implement a stable and powerful multi-tone harmonic balance engine with support to autonomous circuits.

Another goal will be to integrate YalRF with the scikit-rf / openEMS / SignalIntegrity packages.

Example of usage:
```python
import matplotlib.pyplot as plt
from yalrf import Netlist
from yalrf.Analyses import MultiToneHarmonicBalance

net = Netlist('Peltz Oscillator')

# circuit parameters
vcc = 10
r = 200e3
l = 0.5e-3
c = 10e-9
re = 50e3

# VCC
net.add_idc('I1', 'nx', 'gnd', dc=vcc)
net.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

# tank circuit
net.add_resistor('R1', 'nvcc', 'nb', r)
net.add_inductor('L1', 'nvcc', 'nb', l)
net.add_capacitor('C1', 'nvcc', 'nb', c)

# emitter resistance
net.add_resistor('RE', 'ne', 'gnd', re)

# bjts
q1 = net.add_bjt('Q1', 'nb', 'nvcc', 'ne')
q2 = net.add_bjt('Q2', 'nvcc', 'nb', 'ne')

q1.options['Is'] = 1e-16
q1.options['Bf'] = 200
q1.options['Br'] = 1
q2.options = q1.options.copy()

numharmonics = 10
freq = 80e3
V0 = 0.1

hb = MultiToneHarmonicBalance('HB1')
hb.options['maxiter'] = 100

converged, freqs, Vf, _, _ = hb.run_oscillator(net, freq, numharmonics, V0, 'nb')

hb.print_v('nb')
hb.plot_v('nb')
plt.show()
```

# Dependencies

`conda install python numpy scipy matplotlib`

# Acknowledgements

The API for netlist description was based on ahkab. Many device models currently used were referenced from the Qucs documentation.

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
