import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import MultiToneHarmonicBalance, HarmonicBalance
import time

y = YalRF('Diode Testbench')

i1 = y.add_iac('I1', 'nx', 'gnd', ac=1, freq=1.1e6)
g1 = y.add_gyrator('G1', 'nx', 'nz', 'gnd', 'gnd', 1)

i2 = y.add_iac('I2', 'ny', 'gnd', ac=1, freq=0.9e6)
g2 = y.add_gyrator('G2', 'ny', 'n1', 'nz', 'gnd', 1)

r1 = y.add_resistor('R1', 'n1', 'n2', 100)
d1 = y.add_diode('D1', 'n2', 'gnd')

d1.options['Is'] = 1e-15
d1.options['N'] = 1
d1.options['Area'] = 1

begin = time.time()

hb = MultiToneHarmonicBalance('HB1', [1.1e6, 0.9e6], [3, 3])

converged, freqs, Vf, _, _ = hb.run(y)

end = time.time()

hb.print_v('n2')
hb.plot_v('n2')
plt.show()

print('Running time: {}'.format(end-begin))
print('HB problem size: {}'.format(hb.V.shape))
print('Number of frequency bins: {}'.format(freqs.shape))
