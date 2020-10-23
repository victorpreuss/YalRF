import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance

y = YalRF('Diode Testbench')

v1 = y.add_vsource('V1', 'n1', 'gnd', dc=0, ac=1.0)
r1 = y.add_resistor('R1', 'n1', 'n2', 50)
r2 = y.add_resistor('R1', 'n2', 'n3', 50)
# c1 = y.add_capacitor('C1', 'n2', 'gnd', 1e-6)

d1 = y.add_diode('D1', 'n3', 'n4')
d1.options['Is'] = 1e-15
d1.options['N'] = 1
d1.options['Area'] = 1
d1.options['Rs'] = 1e-6

d2 = y.add_diode('D2', 'n4', 'gnd')
d2.options['Is'] = 1e-15
d2.options['N'] = 1
d2.options['Area'] = 1

# remember to always use the harmonics number to a power of 2
hb = HarmonicBalance('HB1', 1e6, 8)

x = hb.run(y)

# print(x)
