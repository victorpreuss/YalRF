import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF

y = YalRF('BJT AC Testbench')

v1 = y.add_vdc('V1', 'nc', 'gnd', 1)

c1 = y.add_capacitor('C1', 'ny', 'nb', 1e-6)
i1 = y.add_iac('I1', 'ny', 'gnd', 1)

l1 = y.add_inductor('L1', 'nx', 'nb', 1e-3)
v2 = y.add_vdc('V2', 'nx', 'gnd', 0.75)

q1 = y.add_bjt('Q1', 'nb', 'nc', 'gnd')
q1.options['Is'] = 8.11e-14
q1.options['Nf'] = 1
q1.options['Nr'] = 1
q1.options['Ikf'] = 0.5
q1.options['Ikr'] = 0.225
q1.options['Vaf'] = 113
q1.options['Var'] = 24
q1.options['Ise'] = 1.06e-11
q1.options['Ne'] = 2
q1.options['Isc'] = 0
q1.options['Nc'] = 2
q1.options['Bf'] = 205
q1.options['Br'] = 4

q1.options['Cje'] = 2.95e-11
q1.options['Cjc'] = 1.52e-11
q1.options['Cjs'] = 0.

dc1 = y.add_dc_analysis('DC1')
xdc = y.run('DC1')

ac1 = y.add_ac_analysis('AC1', start=10e6, stop=10e9, numpts=30, sweeptype='logarithm')
xac = y.run('AC1', xdc)

freqs = y.get_freqs('AC1')
vb = y.get_voltage('AC1', 'nb')
vc = y.get_voltage('AC1', 'nc')

ib = i1.ac
Rin = 1. / np.real(ib / vb)
Cin = np.imag(ib / vb) / (2. * np.pi * freqs)

plt.figure(figsize=(12,4))

plt.subplot(121)
plt.ticklabel_format(useOffset=False)
plt.semilogx(freqs, Rin * 1e-6)
plt.title('BJT AC Testbench')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Rin [MOhms]')
plt.grid()

plt.subplot(122)
plt.semilogx(freqs, Cin * 1e12)
plt.title('BJT AC Testbench')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Cin [pF]')
plt.grid()
plt.show()


