import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import MultiToneHarmonicBalance

y = YalRF('Diode Testbench')

i1 = y.add_iac('I1', 'nx', 'gnd', ac=1, freq=0.9e6)
g1 = y.add_gyrator('G1', 'nx', 'nz', 'gnd', 'gnd', 1)

i2 = y.add_iac('I2', 'ny', 'gnd', ac=1, freq=1.1e6)
g2 = y.add_gyrator('G2', 'ny', 'n1', 'nz', 'gnd', 1)

r1 = y.add_resistor('R1', 'n1', 'n2', 100)

d1 = y.add_diode('D1', 'n2', 'gnd')

d1.options['Is'] = 1e-15
d1.options['N'] = 1
d1.options['Area'] = 1

hb = MultiToneHarmonicBalance('HB1', [0.9e6, 1.1e6], [10, 10])

converged, freqs, Vf, time, Vt = hb.run(y)

#hb.plot_v('n1')
#hb.plot_v('n2')
#plt.show()
