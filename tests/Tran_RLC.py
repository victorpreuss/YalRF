import numpy as np
import matplotlib.pyplot as plt

import setup
from Xyce import getXyceData
from yalrf import YalRF

y = YalRF('RLC Transient')

y.add_vpulse('V1', 'n1', 'gnd', v1=0, v2=1, tstart=0, tstop=1)
y.add_resistor('R1', 'n1', 'n2', 50)
y.add_capacitor('C1', 'n2', 'gnd', 1.074e-12)
y.add_inductor('L1', 'n2', 'n3', 7.723e-9)
y.add_capacitor('C2', 'n3', 'gnd', 7.014e-12)
y.add_resistor('R2', 'n3', 'gnd', 50)

tr1 = y.add_tran_analysis('TR1', 2e-9, maxtstep=5e-12)
#tr1.options['max_iterations'] = 500
#tr1.options['reltol'] = 1e-6

y.run('TR1')

t = y.get_time('TR1')
vn2 = y.get_voltage('TR1', 'n2')
vn3 = y.get_voltage('TR1', 'n3')
ir1 = y.get_itran('TR1', 'R1')
il1 = y.get_itran('TR1', 'L1')
ic1 = y.get_itran('TR1', 'C1')

plt.figure(figsize=(12,4))

plt.subplot(121)
plt.plot(t, ir1 * 1e3)
plt.plot(t, il1 * 1e3)
plt.plot(t, ic1 * 1e3)
plt.xlabel('Time [s]')
plt.ylabel('Current [mA]')
plt.legend(['R1', 'L1', 'C1'])
plt.grid()

plt.subplot(122)
plt.plot(t, vn2)
plt.plot(t, vn3)
plt.title('RLC Transient')
plt.xlabel('Time [s]')
plt.ylabel('Voltage [V]')
plt.legend(['vn2', 'vn3'])
plt.grid()
plt.show()