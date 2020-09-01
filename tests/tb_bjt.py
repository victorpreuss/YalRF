import numpy as np
import matplotlib.pyplot as plt

import setup
from yarf import Yarf

y = Yarf('BJT AC Testbench')

v1 = y.add_vdc('V1', 'nc', 'gnd', 1)

c1 = y.add_capacitor('C1', 'ny', 'nb', 1e-6)
i1 = y.add_isource('I1', 'ny', 'gnd', dc=0, ac=100e-6)

l1 = y.add_inductor('L1', 'nx', 'nb', 1e-6)
v2 = y.add_vdc('V2', 'nx', 'gnd', 0.6)

q1 = y.add_bjt('Q1', 'nb', 'nc', 'gnd')
q1.options['Cje'] = 1e-12
q1.options['Cjc'] = 1e-12
q1.options['Cjs'] = 1e-12

dc1 = y.add_dc_analysis('DC1')
xdc = y.run('DC1')

ac1 = y.add_ac_analysis('AC1', start=10e6, stop=10e9, numpts=30, sweeptype='logarithm')
xac = y.run('AC1', xdc)

freqs = y.get_freqs('AC1')
vb = y.get_voltage('AC1', 'nb')
vc = y.get_voltage('AC1', 'nc')

Cin = np.imag(100e-6 / vb) / (2. * np.pi * freqs)

print('{}'.format(Cin[0]))
print('{}'.format(Cin[-1]))

plt.figure(1)
plt.semilogx(freqs, Cin)
plt.grid()
plt.title('AC Simulation')
plt.xlabel('Frequencies')
plt.ylabel('Cin')
plt.show()


