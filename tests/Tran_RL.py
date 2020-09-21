import numpy as np
import matplotlib.pyplot as plt

import setup
from Xyce import getXyceData
from yalrf import YalRF

y = YalRF('RL Transient')

y.add_vpulse('V1', 'n1', 'gnd', v1=0, v2=1, tstart=0, tstop=10e-9)
y.add_resistor('R1', 'n1', 'n2', 1e3)
y.add_inductor('L1', 'n2', 'gnd', 1e-6)

tr1 = y.add_tran_analysis('TR1', tstop=20e-9, maxtstep=200e-12)
#tr1.options['max_iterations'] = 500
#tr1.options['reltol'] = 1e-6

y.run('TR1')

t = y.get_time('TR1')
i = y.get_itran('TR1', 'R1')

plt.figure(figsize=(10,5))
plt.plot(t, i * 1e3)
plt.title('RL Transient')
plt.xlabel('Time [s]')
plt.ylabel('I [mA]')
plt.grid()
plt.show()
