import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("../yarf")

from yarf import Yarf

y = Yarf("Hello World!")

y.add_vsource('V1', 'n1', 'gnd', dc=0, ac=1)
y.add_resistor('R1', 'n1', 'n2', 50.)
y.add_inductor('L1', 'n2', 'n3', 1e-6)
y.add_vdc('V2', 'n3', 'gnd', 0.75)
d1 = y.add_diode('D1', 'n2', 'gnd')
d1.options['Cj0'] = 1e-12

ac = y.add_ac_analysis('AC1', start=10e6, stop=10e9, numpts=300, sweeptype='logarithm')

sol = y.run('AC1')

freq = y.get_freqs('AC1')
vout = y.get_voltage('AC1', 'n2')

print('{}'.format(np.abs(vout[0])))
print('{}'.format(np.abs(vout[-1])))

plt.figure(1)
plt.semilogx(freq, np.abs(vout))
plt.grid()
plt.title('AC Simulation')
plt.xlabel('Frequencies')
plt.ylabel('Vout')
plt.show()
