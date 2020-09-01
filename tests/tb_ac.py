import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../yarf")

from yarf import Yarf

y = Yarf("Hello World!")

y.add_vsource('V1', 'n1', 'gnd', dc=0, ac=1)

y.add_resistor('R1', 'n1', 'n2', 50.)
y.add_inductor('L1', 'n2', 'n3', 0.245894)
y.add_capacitor('C1', 'n3', 'n4', 1.03013e-07)
y.add_inductor('L2', 'n4', 'gnd', 9.83652e-05)
y.add_capacitor('C2', 'n4', 'gnd', 0.000257513)
y.add_inductor('L3', 'n4', 'n5', 0.795775)
y.add_capacitor('C3', 'n5', 'n6', 3.1831e-08)
y.add_inductor('L4', 'n6', 'gnd', 9.83652e-05)
y.add_capacitor('C4', 'n6', 'gnd', 0.000257513)
y.add_capacitor('C5', 'n7', 'n8', 1.03013e-07)
y.add_inductor('L5', 'n6', 'n7', 0.245894)
y.add_resistor('R2', 'n8', 'gnd', 50.)

ac = y.add_ac_analysis('AC1', start=970, stop=1030, numpts=50)

sol = y.run('AC1')

freq = y.get_freqs('AC1')
vout = y.get_voltage('AC1', 'n8')

plt.figure(1)
plt.plot(freq, np.abs(vout))
plt.grid()
plt.title('AC Simulation')
plt.xlabel('Frequencies')
plt.ylabel('Vout')
plt.show()
