import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance, MultiToneHarmonicBalance

y = Netlist('Differential Amplifier')

# circuit parameters
ibias = 50e-3
vbias = 2.0
vcc = 5
vin = 60e-3
rload = 50

# VCC
i1 = y.add_idc('I1', 'nx', 'gnd', dc=vcc)
g1 = y.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

# vin1
i2 = y.add_iac('I2', 'ny', 'gnd', ac=vin, phase=-90, freq=10e6)
i3 = y.add_idc('I3', 'ny', 'gnd', dc=vbias)
g2 = y.add_gyrator('G2', 'ny', 'nb1', 'gnd', 'gnd', 1)

# vin2
i4 = y.add_iac('I4', 'nz', 'gnd', ac=vin, phase=+90, freq=10e6)
i5 = y.add_idc('I5', 'nz', 'gnd', dc=vbias)
g3 = y.add_gyrator('G3', 'nz', 'nb2', 'gnd', 'gnd', 1)

# collector loads
r1 = y.add_resistor('R1', 'nvcc', 'nc1', rload)
r2 = y.add_resistor('R2', 'nvcc', 'nc2', rload)

# differential pair
q1 = y.add_bjt('Q1', 'nb1', 'nc1', 'ne')
q2 = y.add_bjt('Q2', 'nb2', 'nc2', 'ne')

q3 = y.add_bjt('Q3', 'nbx', 'ne', 'gnd')
q4 = y.add_bjt('Q4', 'nbx', 'nbx', 'gnd')

i6 = y.add_idc('I6', 'nbx', 'gnd', dc=ibias)

q1.options['Is'] = 1e-15
q1.options['Bf'] = 100
q1.options['Vaf'] = 20

q2.options = q1.options
q3.options = q1.options
q4.options = q1.options

# run harmonic balance
# hb = HarmonicBalance('HB1', 10e6, 10)
hb = MultiToneHarmonicBalance('HB1', 10e6, 20)
hb.options['maxiter'] = 100
converged, freqs, Vf, time, Vt = hb.run(y)

hb.print_v('nc1')
hb.plot_v('nc1')
# hb.plot_v('nc2')
plt.show()
