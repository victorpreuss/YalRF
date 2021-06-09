import time
import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import MultiToneHarmonicBalance


y = Netlist('Van Der Pol Oscillator')

# circuit parameters
C = 1
L = 1
alpha = 2.5

y.add_inductor('L', 'np', 'gnd', L)
y.add_capacitor('C', 'np', 'gnd', C)
r = y.add_resistor('R', 'np', 'gnd', -1/alpha)

# this element is i = alpha * v^3
x = y.add_cubicnl('X1', 'np', 'gnd', alpha)

numharmonics = 20
freq = 0.13
V0 = 1.2

hb = MultiToneHarmonicBalance('HB1')
hb.options['maxiter'] = 100

# converged, freqs, Vf, _, _ = hb.run_oscillator(y, freq, numharmonics, V0, 'np')

# calculate the inductor current as: I = V / (jwL)
# v = hb.get_v('np')
# i = np.zeros(v.shape, dtype=complex)
# for idx in range(1,len(freqs)):
#     i[idx] = v[idx] / (1j * (2 * np.pi * freqs[idx]) * L)

# t, vt = hb.convert_to_time(v)
# t, it = hb.convert_to_time(i)

# hb.print_v('np')
# hb.plot_v('np')

# plt.figure()
# plt.plot(t, it)

# plt.figure()
# plt.plot(t, alpha * (vt**3))

plt.figure(figsize=(8,6))
plt.xlabel('Voltage [V]')
plt.ylabel('Current [A]')
plt.xlim((-1.5,1.5))
plt.ylim((-2,2))
plt.grid()

for alpha in [0.1, 1.0, 1.8, 2.5]:

    r.R = -1/alpha
    x.alpha = alpha

    converged, freqs, Vf, _, _ = hb.run_oscillator(y, freq, numharmonics, V0, 'np')

    # calculate the inductor current as: I = V / (jwL)
    v = hb.get_v('np')
    i = np.zeros(v.shape, dtype=complex)
    for idx in range(1,len(freqs)):
        i[idx] = v[idx] / (1j * (2 * np.pi * freqs[idx]) * L)

    t, vt = hb.convert_to_time(v)
    t, it = hb.convert_to_time(i)

    plt.plot(vt, it, label=r'$\alpha$ = {:.1f}'.format(alpha))

    idx = 0
    plt.arrow(vt[idx], it[idx], vt[idx+1]-vt[idx], it[idx+1]-it[idx], shape='full', lw=0, length_includes_head=True, head_width=.09, facecolor='red')

    idx = int(len(vt) / 8.)
    plt.arrow(vt[idx], it[idx], vt[idx+1]-vt[idx], it[idx+1]-it[idx], shape='full', lw=0, length_includes_head=True, head_width=.09, facecolor='red')

    idx = int(len(vt) / 4.)
    plt.arrow(vt[idx], it[idx], vt[idx+1]-vt[idx], it[idx+1]-it[idx], shape='full', lw=0, length_includes_head=True, head_width=.09, facecolor='red')

    idx = int(3. * len(vt) / 8.)
    plt.arrow(vt[idx], it[idx], vt[idx+1]-vt[idx], it[idx+1]-it[idx], shape='full', lw=0, length_includes_head=True, head_width=.09, facecolor='red')

plt.legend(loc='upper right')

plt.show()


