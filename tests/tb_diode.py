import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../yarf")

from yarf import Yarf

y = Yarf('Diode Testbench')

i1 = y.add_idc('I1', 'n1', 'gnd', 1e-3)
y.add_resistor('R1', 'n1', 'n2', 50)

d1 = y.add_diode('D1', 'n2', 'gnd')
d1.options['Is'] = 1e-15
d1.options['N']  = 1

dc1 = y.add_dc_analysis('DC1')
dc1.options['max_iterations'] = 50

y.run('DC1')

y.print_dc_voltages('DC1')
y.print_dc_currents('DC1')

vd  = []
idc = np.linspace(1e-3, 100e-3, 100)
for i in idc:
    i1.dc = i
    x = y.run('DC1')
    vd.append(x[1,0])

plt.figure(1)
plt.plot(vd, idc)
plt.grid()
plt.title('DC IxV curve of a Diode')
plt.xlabel('Vd [V]')
plt.ylabel('Id [A]')
plt.show()


