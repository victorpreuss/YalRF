import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance

y = YalRF('Full-Wave Rectifier Testbench')

# i1 = y.add_iac('I1', 'vinp', 'gnd', ac=15e-3, phase=+90)
# i2 = y.add_idc('I2', 'vinp', 'gnd', dc=0e-3)

# r1 = y.add_resistor('R1', 'vinp', 'gnd', 100)
# r2 = y.add_resistor('R2', 'voutp', 'voutn', 100)

# d1 = y.add_diode('D1', 'vinp', 'voutp')
# d2 = y.add_diode('D2', 'gnd', 'voutp')
# d3 = y.add_diode('D3', 'voutn', 'vinp')
# d4 = y.add_diode('D4', 'voutn', 'gnd')

i1 = y.add_iac('I1', 'vinp', 'vinn', ac=15e-3, phase=+90)
i2 = y.add_idc('I2', 'vinp', 'vinn', dc=5e-3)

r1 = y.add_resistor('R1', 'vinp', 'vinn', 100)
r2 = y.add_resistor('R2', 'voutp', 'gnd', 100)

d1 = y.add_diode('D1', 'vinp', 'voutp')
d2 = y.add_diode('D2', 'vinn', 'voutp')
d3 = y.add_diode('D3', 'gnd', 'vinp')
d4 = y.add_diode('D4', 'gnd', 'vinn')

d1.options['Is'] = 1e-15
d1.options['N'] = 1
d1.options['Area'] = 1

d2.options['Is'] = 1e-15
d2.options['N'] = 1
d2.options['Area'] = 1

d3.options['Is'] = 1e-15
d3.options['N'] = 1
d3.options['Area'] = 1

d4.options['Is'] = 1e-15
d4.options['N'] = 1
d4.options['Area'] = 1

hb = HarmonicBalance('HB1', 1e6, 10)

V = hb.run(y)

n1 = 2 # voutp
n2 = 0 # voutn

K = hb.numharmonics
freqs = hb.freq * np.linspace(0, K, K+1)
Kk = 2 * (K + 1)

# assemble complex array of spectra for nodes 'voutp' and 'voutn'
vf = np.zeros(K+1, dtype=complex)
for k in range(K+1):
    vf[k] = (V[Kk*n1+2*k+0] + 1j * V[Kk*n1+2*k+1])# - (V[Kk*n2+2*k+0] + 1j * V[Kk*n2+2*k+1])

# compute inverse fourier transform of voltage waveform
S = 8 * K
vt = np.zeros(S)
for s in range(S):
    vt[s] = vf[0].real
    for k in range(1, K+1):
        vt[s] = vt[s] + 2 * (vf[k].real * np.cos(2. * np.pi * k * s / S) -
                             vf[k].imag * np.sin(2. * np.pi * k * s / S))

vt_plot1 = vt.copy()
vf_plot1 = vf.copy()

plt.figure()
plt.subplot(211)
plt.plot(vt_plot1)
plt.grid()
plt.subplot(212)
plt.stem(freqs, abs(vf_plot1), use_line_collection=True, markerfmt='^')
for f, v in zip(freqs, vf_plot1):
    label = "({:.3f}, {:.1f})".format(abs(v), np.degrees(np.angle(v)))
    plt.annotate(label, (f,abs(v)), textcoords="offset points", xytext=(0,10), ha='center')
plt.grid()
plt.show()