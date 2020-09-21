import numpy as np
import matplotlib.pyplot as plt

import setup
from Xyce import getXyceData
from yalrf import YalRF

y = YalRF('RC Transient')

y.add_vpulse('V1', 'n1', 'gnd', v1=0, v2=1, tstart=0, tstop=25e-3, trise=1e-9, tfall=1e-9)
y.add_resistor('R1', 'n1', 'n2', 1e3)
y.add_capacitor('C1', 'n2', 'gnd', 1e-6)

tr1 = y.add_tran_analysis('TR1', tstop=50e-3, maxtstep=500e-6)

y.run('TR1')

t = y.get_time('TR1')
vn2 = y.get_voltage('TR1', 'n2')

plt.figure(figsize=(10,5))
plt.plot(t, vn2)
plt.title('RC Transient')
plt.xlabel('Time [s]')
plt.ylabel('Vc [V]')
plt.grid()
plt.show()