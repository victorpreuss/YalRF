import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF

y = YalRF('Diode Testbench')

v1 = y.add_vdc('V1', 'n1', 'gnd', 0.5)

d1 = y.add_diode('D1', 'n1', 'gnd')
d1.options['Is'] = 1e-15
d1.options['N'] = 1
d1.options['Rs'] = 0.1
d1.options['Area'] = 2

dc1 = y.add_dc_analysis('DC1')
#dc1.options['reltol'] = 1e-3

vdc = np.linspace(0.5, 1, 100)
Id  = []
for v in vdc:
    v1.dc = v
    x = y.run('DC1')
    Id.append(-x[1,0])

plt.figure()
plt.plot(vdc, Id)
plt.grid()
plt.title('DC IxV curve of a Diode')
plt.xlabel('Vd [V]')
plt.ylabel('Id [A]')
plt.show()


