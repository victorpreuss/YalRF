import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF

y = YalRF('RC Transient')

y.add_vstep('V1', 'n1', 'gnd', 1.0, tstart=0, tstop=1)
y.add_resistor('R1', 'n1', 'n2', 50)
y.add_capacitor('C1', 'n2', 'gnd', 1.074e-12)
y.add_inductor('L1', 'n2', 'n3', 7.723e-9)
c = y.add_capacitor('C2', 'n3', 'gnd', 7.014e-12)
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

plt.figure(figsize=(10,5))
plt.plot(t, ir1)
plt.plot(t, il1)
plt.plot(t, ic1)
plt.legend(['R1', 'L1', 'C1'])
plt.grid()

plt.figure(figsize=(10,5))
plt.plot(t, vn2)
plt.plot(t, vn3)
plt.title('RC Transient')
plt.xlabel('Time [s]')
plt.ylabel('Vc [V]')
plt.legend(['vn2', 'vn3'])
plt.grid()
plt.show()


