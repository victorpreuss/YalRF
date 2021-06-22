import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance, MultiToneHarmonicBalance


y = Netlist('Oscillator')

# circuit parameters
vcc = 10
r = 200e3
l = 0.5e-3
c = 10e-9
re = 50e3
freq = 0

# VCC
y.add_idc('I1', 'nx', 'gnd', dc=vcc)
y.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

# tank circuit
y.add_resistor('R1', 'nvcc', 'nb', r)
y.add_inductor('L1', 'nvcc', 'nb', l)
y.add_capacitor('C1', 'nvcc', 'nb', c)

# emitter resistance
y.add_resistor('RE', 'ne', 'gnd', re)

# bjts
q1 = y.add_bjt('Q1', 'nb', 'nvcc', 'ne')
q2 = y.add_bjt('Q2', 'nvcc', 'nb', 'ne')

q1.options['Is'] = 1e-16
q1.options['Bf'] = 200
q1.options['Br'] = 1
q2.options = q1.options.copy()

numharmonics = 10
freq = 80e3
V0 = 0.1

hb = MultiToneHarmonicBalance('HB1')
hb.options['maxiter'] = 100

converged, freqs, Vf, _, _ = hb.run_oscillator(y, freq, numharmonics, V0, 'nb')

hb.print_v('nb')
hb.plot_v('nb')
# hb.plot_v('ne')
plt.show()


