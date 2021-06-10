import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import MultiToneHarmonicBalance

y = Netlist('Oscillator')

# circuit parameters
vcc = 5
c1 = 5e-9
c2 = 5e-9
c3 = 5e-9
r1 = 15e3
r2 = 15e3
r3 = 2.2e3
r4 = 5e3
re = 5

# VCC
y.add_idc('I1', 'nx', 'gnd', dc=vcc)
y.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

# bias resistors
y.add_resistor('R4', 'nvcc', 'nc', r4)
y.add_resistor('RE', 'ne', 'gnd', re)

# RC feedback network
y.add_capacitor('C1', 'nc', 'n2', c1)
y.add_capacitor('C2', 'n2', 'nb', c2)
y.add_capacitor('C3', 'n1', 'gnd', c3)

y.add_resistor('R1', 'nc', 'n1', r1)
y.add_resistor('R2', 'n1', 'nb', r2)
y.add_resistor('R3', 'n2', 'gnd', r3)

# bjt
q1 = y.add_bjt('Q1', 'nb', 'nc', 'ne')

q1.options['Is'] = 10.2e-15
q1.options['Bf'] = 301
q1.options['Br'] = 4
q1.options['Ne'] = 2
q1.options['Nc'] = 2
q1.options['Ise'] = 5.82e-12
q1.options['Vaf'] = 121
q1.options['Var'] = 24
q1.options['Ikf'] = 60.7e-3
q1.options['Ikr'] = 0.15

q1.options['Cje'] = 26.8e-12
q1.options['Vje'] = 1.1
q1.options['Mje'] = 0.5
q1.options['Cjc'] = 8.67e-12
q1.options['Vjc'] = 0.3
q1.options['Mjc'] = 0.3
q1.options['Xcjc'] = 1
q1.options['Cjs'] = 0
q1.options['Vjs'] = 0.75
q1.options['Mjs'] = 0
q1.options['Fc'] = 0.5
q1.options['Tf'] = 427e-12
q1.options['Tr'] = 50.3e-9

numharmonics = 10
freq = 3.8e3
V0 = 0.5

hb = MultiToneHarmonicBalance('HB1')
hb.options['maxiter'] = 100

converged, freqs, Vf, _, _ = hb.run_oscillator(y, freq, numharmonics, V0, 'nc')

hb.plot_v('nc')
plt.show()

