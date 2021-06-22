import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import AC, MultiToneHarmonicBalance

# y = Netlist('Oscillator')

# # circuit parameters
# vcc = 5
# vbias = 1
# ibias = 1e-3
# rl = 850
# fosc = 50e6
# l = 50e-9
# ct = 1 / (l * (2 * np.pi * fosc) ** 2)

# n = 0.2 # capactive divider ratio (feedback)
# c1 = ct / (1 - n)
# c2 = ct / n

# rs = l ** 2 * (2 * np.pi * fosc) ** 2 / rl

# numharmonics = 20

# f0 = 1 / (2 * np.pi * np.sqrt(l * ct))
# V0 = 2 * ibias * rl * (1 - n)

# # VCC
# y.add_idc('I1', 'nx1', 'gnd', dc=vcc)
# y.add_gyrator('G1', 'nx1', 'nvcc', 'gnd', 'gnd', 1)

# # Vbias
# y.add_idc('I2', 'ny1', 'gnd', dc=vbias)
# y.add_gyrator('G2', 'ny1', 'nb', 'gnd', 'gnd', 1)

# # Ibias
# y.add_idc('I3', 'gnd', 'ne', dc=ibias)

# # passives
# # y.add_resistor('Rl', 'nvcc', 'nc', rl)
# y.add_resistor('Rs', 'nvcc', 'nvccx', rs)
# y.add_inductor('L1', 'nvccx', 'nc', l)
# C1 = y.add_capacitor('C1', 'nc', 'ne', c1)
# C2 = y.add_capacitor('C2', 'ne', 'gnd', c2)

# # bjts
# q1 = y.add_bjt('Q1', 'nb', 'nc', 'ne')

# q1.options['Is'] = 1e-15
# q1.options['Bf'] = 100

####################################################

vin = 10e-3
fin = 50e6

vcc = 5
vbias = 1
ibias = 1e-3

net2 = Netlist('Large Gm')

# VCC
net2.add_idc('I1', 'nx', 'gnd', dc=vcc)
net2.add_gyrator('G1', 'nx', 'nc', 'gnd', 'gnd', 1)

net2.add_idc('I2', 'ny', 'gnd', dc=vbias)
iin = net2.add_iac('I3', 'ny', 'gnd', ac=vin, freq=fin)
net2.add_gyrator('G2', 'ny', 'nb', 'gnd', 'gnd', 1)

net2.add_idc('I4', 'gnd', 'ne', dc=ibias)
net2.add_capacitor('C1', 'ne', 'gnd', 1e-6)

q1 = net2.add_bjt('Q1', 'nb', 'nc', 'ne')

q1.options['Is'] = 1e-15
q1.options['Bf'] = 100

hb = MultiToneHarmonicBalance('HB1', fin, 20)
hb.options['maxiter'] = 100

vi = np.arange(100e-6, 251e-3, 5e-3)
ic = np.zeros(vi.shape)
V0 = None
i = 0
for vin in vi:

    # update input voltage
    iin.ac = vin

    # run harmonic balance
    converged, freqs, Vf, _, _ = hb.run(net2, V0)
    V0 = hb.V

    # get time-domain waveform from BJT and convert to frequency
    ic_td = np.array(q1.Ic, dtype=complex)
    ic_fd = hb.DFT @ ic_td
    ic[i] = np.abs(ic_fd[1] + 1j * ic_fd[2])

    i += 1

from scipy.constants import k, e

gm = ibias / (k * 300 / e)

plt.figure()
plt.plot(vi * 1e3, ic / vi / gm, label='Simulation')
plt.plot(vi[vi>=50e-3] * 1e3, 2 * ibias / vi[vi>=50e-3] / gm, label='Approximation')
plt.grid()
plt.legend()
plt.xlabel('$Vin$ [mV]')
plt.ylabel('$G_m / g_m$')
plt.show()
