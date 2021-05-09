import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance, MultiToneHarmonicBalance

y = YalRF('Differential Amplifier')

# circuit parameters
ibias = 10e-3
vbias = 1.5
vcc = 5
vin = 60e-3

# VCC
i1 = y.add_idc('I1', 'nx', 'gnd', dc=vcc)
g1 = y.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

# vin1
i2 = y.add_iac('I2', 'ny', 'gnd', ac=vin, phase=0, freq=10e6)
i3 = y.add_idc('I3', 'ny', 'gnd', dc=vbias)
g2 = y.add_gyrator('G2', 'ny', 'nb1', 'gnd', 'gnd', 1)

# vin2
i4 = y.add_iac('I4', 'nz', 'gnd', ac=vin, phase=+180, freq=10e6)
i5 = y.add_idc('I5', 'nz', 'gnd', dc=vbias)
g3 = y.add_gyrator('G3', 'nz', 'nb2', 'gnd', 'gnd', 1)

# differential pair
q1 = y.add_bjt('Q1', 'nb1', 'nc1', 'ne')
q2 = y.add_bjt('Q2', 'nb2', 'nc2', 'ne')

# active loads
q5 = y.add_bjt('Q5', 'nc1', 'nc1', 'nvcc', ispnp=True)
q6 = y.add_bjt('Q6', 'nc1', 'nc2', 'nvcc', ispnp=True)

# current mirror
q3 = y.add_bjt('Q3', 'nbx', 'ne', 'gnd')
q4 = y.add_bjt('Q4', 'nbx', 'nbx', 'gnd')

i6 = y.add_idc('I6', 'nbx', 'gnd', dc=ibias)

q1.options['Is'] = 1e-15
q1.options['Bf'] = 100
q1.options['Vaf'] = 60

q2.options = q1.options
q3.options = q1.options
q4.options = q1.options

# pnp transistors
q5.options = q1.options
q6.options = q1.options

dc1 = y.add_dc_analysis('DC1')
xdc = y.run('DC1')
y.print_dc_voltages('DC1')

# ac1 = y.add_ac_analysis('AC1', start=10e3, stop=100e6, numpts=1000, sweeptype='linear')
# xac = y.run('AC1', xdc)


