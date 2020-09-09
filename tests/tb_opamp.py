import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF

y = YalRF('Simple inverting amplifier')

y.add_vsine('V1', 'n1', 'gnd', dc=1, ac=1, freq=1e9, phase=0)
y.add_opamp('OP1', 'gnd', 'n2', 'n3', G=100e3, Vmax=100)
y.add_resistor('R1', 'n1', 'n2', 50)
y.add_resistor('R2', 'n2', 'n3', 100)
y.add_capacitor('C1', 'n2', 'n3', 1e-12)

y.add_tran_analysis('TR1', 2e-9, maxtstep=5e-12)
y.add_ac_analysis('AC1', start=10e6, stop=10e9, numpts=100)

y.run('TR1')
y.run('AC1')

# get transient data
t = y.get_time('TR1')
vn1 = y.get_voltage('TR1', 'n1')
vn2 = y.get_voltage('TR1', 'n2')
vn3 = y.get_voltage('TR1', 'n3')
ir1 = y.get_itran('TR1', 'R1')

# get ac data
f = y.get_freqs('AC1')
G = 20. * np.log10(np.abs(y.get_voltage('AC1', 'n3') / y.get_voltage('AC1', 'n1')))

plt.figure()
plt.semilogx(f / 1e6, G)
plt.title('Inverting Amplifier Gain')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Gain [dB]')
plt.grid()

plt.figure(figsize=(12,5))

plt.subplot(121)
plt.plot(t, vn1)
plt.plot(t, vn2)
plt.plot(t, vn3)
plt.title('Voltages')
plt.xlabel('Time [s]')
plt.ylabel('Voltage [V]')
plt.legend(['vn1', 'vn2', 'vn3'])
plt.grid()

plt.subplot(122)
plt.plot(t, ir1 * 1e3)
plt.title('Currents')
plt.xlabel('Time [s]')
plt.ylabel('Current [mA]')
plt.legend(['IR1'])
plt.grid()

plt.show()