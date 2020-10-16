import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF
from yalrf.Analyses import HarmonicBalance

y = YalRF('Diode Testbench')

v1 = y.add_vac('V1', 'n1', 'gnd', 0.5)
r1 = y.add_resistor('R1', 'n1', 'n2', 50)

d1 = y.add_diode('D1', 'n2', 'n3')
d1.options['Is'] = 1e-15
d1.options['N'] = 1
d1.options['Area'] = 1

d2 = y.add_diode('D2', 'n3', 'gnd')
d2.options['Is'] = 1e-15
d2.options['N'] = 1
d2.options['Area'] = 1

hb = HarmonicBalance('HB1', 1e6, 3)

x = hb.run(y)

# print(x)
