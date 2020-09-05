import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import v2_yalrf

y = v2_yalrf('Diode Testbench')

i1 = y.add_vdc('I1', 'n1', 'gnd', 0.3)
#y.add_resistor('R1', 'n1', 'n2', 1e3)

d1 = y.add_diode('D1', 'n1', 'gnd')
d1.options['Is'] = 1e-15
d1.options['N']  = 1

dc1 = y.add_dc_analysis('DC1')
#dc1.options['reltol'] = 1e-9
dc1.options['iabstol'] = 1e-12

vd  = []
idc = np.linspace(-10, 1, 100)
for i in idc:
    i1.dc = i
    x = y.run('DC1')
    vd.append(-x[1,0])

plt.figure(1)
#plt.plot(vd, idc)
plt.plot(idc, vd)
plt.grid()
plt.title('DC IxV curve of a Diode')
plt.xlabel('Vd [V]')
plt.ylabel('Id [A]')
plt.show()


