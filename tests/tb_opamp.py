import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF

y = YalRF('Simple inverting amplifier')

y.add_vsine('V1', 'n1', 'gnd', dc=1, ac=1, freq=1e9, phase=0)
y.add_opamp('OP1', 'gnd', 'n2', 'n3', G=100e3, Vmax=100)
y.add_resistor('R1', 'n1', 'n2', 50)
y.add_resistor('R2', 'n2', 'n3', 150)

tr1 = y.add_tran_analysis('TR1', 2e-9, maxtstep=5e-12)
#tr1.options['max_iterations'] = 500
#tr1.options['reltol'] = 1e-6

y.run('TR1')

t = y.get_time('TR1')
vn1 = y.get_voltage('TR1', 'n1')
vn2 = y.get_voltage('TR1', 'n2')
vn3 = y.get_voltage('TR1', 'n3')
ir1 = y.get_itran('TR1', 'R1')

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