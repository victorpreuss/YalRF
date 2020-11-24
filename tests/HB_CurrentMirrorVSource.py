import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance

y = YalRF('Differential Amplifier')

vin = 10e-3

# VCC
i1 = y.add_idc('I1', 'nx', 'gnd', dc=5)
g1 = y.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

# vin1
i2 = y.add_iac('I2', 'ny', 'gnd', ac=vin, phase=0)
i3 = y.add_idc('I3', 'ny', 'gnd', dc=2.5)
g2 = y.add_gyrator('G2', 'ny', 'nb1', 'gnd', 'gnd', 1)

# vin2
i4 = y.add_iac('I4', 'nz', 'gnd', ac=vin, phase=180)
i5 = y.add_idc('I5', 'nz', 'gnd', dc=2.5)
g3 = y.add_gyrator('G3', 'nz', 'nb2', 'gnd', 'gnd', 1)

# collector loads
r1 = y.add_resistor('R1', 'nvcc', 'nc1', 100)
r2 = y.add_resistor('R2', 'nvcc', 'nc2', 100)

# differential pair
q1 = y.add_bjt('Q1', 'nb1', 'nc1', 'ne')
q2 = y.add_bjt('Q2', 'nb2', 'nc2', 'ne')

# emitter resistor
# r3 = y.add_resistor('R3', 'ne', 'gnd', 100)
i6 = y.add_idc('I6', 'gnd', 'ne', dc=50e-3)

q1.options['Is'] = 1e-15
q1.options['Bf'] = 100
q1.options['Vaf'] = 20

hb = HarmonicBalance('HB1', 10e6, 5)

freqs, Vf, Vt = hb.run(y)

n = 7
plt.figure()
plt.subplot(211)
plt.plot(Vt[n,:])
plt.plot(Vt[n-1,:])
plt.grid()
plt.subplot(212)
plt.stem(freqs, abs(Vf[n,:]), use_line_collection=True, markerfmt='^')
for f, v in zip(freqs, Vf[n,:]):
    label = "({:.3f}, {:.1f})".format(abs(v), np.degrees(np.angle(v)))
    plt.annotate(label, (f,abs(v)), textcoords="offset points", xytext=(0,10), ha='center')
plt.grid()
plt.show()
