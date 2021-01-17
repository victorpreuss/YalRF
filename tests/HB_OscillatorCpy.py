import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance

y = Netlist('Oscillator')

# circuit parameters
vcc = 5
c1 = 3.3e-12
c2 = 3.9e-12
r1 = 3e3
r2 = 6.8e3
re = 510
l1 = 10e-9

vin = 2.797
freq = 1.205e9

# VCC
y.add_idc('I1', 'nx', 'gnd', dc=vcc)
y.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

# vin oscprobe
Voscprobe = y.add_iac('I2', 'ny', 'gnd', ac=vin, phase=0)
y.add_gyrator('G2', 'ny', 'nz', 'gnd', 'gnd', 1)

# ideal harmonic filter of oscprobe
Zoscprobe = y.add_idealharmonicfilter('X1', 'nz', 'nind', freq)

# bias resistors
y.add_resistor('R1', 'nvcc', 'nb', r1)
y.add_resistor('R2', 'nb', 'gnd', r2)

# emitter resistance
y.add_resistor('RE', 'ne', 'gnd', re)

# capacitor feedback network
y.add_capacitor('C1', 'nb', 'ne', c1)
y.add_capacitor('C2', 'ne', 'gnd', c2)

# resonating inductor
y.add_inductor('L1', 'nind', 'gnd', l1)

# dc feed and dc block (TODO: improve)
y.add_inductor('Lfeed', 'nvcc', 'nc', 1e-3)
y.add_capacitor('Cblk1', 'nb', 'nind', 1e-6)

# load connection
y.add_capacitor('Cblk2', 'nc', 'nl', 1e-6)
y.add_resistor('RL', 'nl', 'gnd', 50)

# bjt
q1 = y.add_bjt('Q1', 'nb', 'nc', 'ne')

q1.options['Is'] = 1e-15
q1.options['Bf'] = 100
q1.options['Br'] = 5
q1.options['Vaf'] = 20
q1.options['Var'] = 10

hb = HarmonicBalance('HB1', freq, 7)
hb.options['maxiter'] = 50

converged, freqs, Vf, time, Vt = hb.run(y)

hb.plot_v('nb')
hb.plot_v('nl')
plt.show()


