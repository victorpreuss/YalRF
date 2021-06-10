import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance, MultiToneHarmonicBalance


y = Netlist('Oscillator')

# circuit parameters
vcc = 5
c1 = 3.3e-12
c2 = 3.9e-12
r1 = 3e3
r2 = 6.8e3
re = 1.5e3
l1 = 6.944e-9

# VCC
y.add_idc('I1', 'nx', 'gnd', dc=vcc)
y.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

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
q1.options['Vaf'] = 60
q1.options['Var'] = 20

numharmonics = 30
freq = 1.2e9
V0 = 1

hb = MultiToneHarmonicBalance('HB1')
hb.options['maxiter'] = 100

converged, freqs, Vf, _, _ = hb.run_oscillator(y, freq, numharmonics, V0, 'nind')

hb.print_v('nind')
hb.plot_v('nind')
# hb.plot_v('nl')
plt.show()


