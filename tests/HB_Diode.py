import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance2

y = YalRF('Diode Testbench')

v1 = y.add_vsource('V1', 'n1', 'gnd', dc=0, ac=3, phase=-90)
r1 = y.add_resistor('R1', 'n2', 'gnd', 100)
#c1 = y.add_capacitor('C1', 'n3', 'gnd', 1e-6)

d1 = y.add_diode('D1', 'n1', 'n2')
d1.options['Is'] = 1e-15
d1.options['N'] = 1
d1.options['Area'] = 1

# d2 = y.add_diode('D2', 'n1', 'n3')
# d2.options['Is'] = 1e-15
# d2.options['N'] = 1
# d2.options['Area'] = 1

# d3 = y.add_diode('D3', 'gnd', 'n2')
# d3.options['Is'] = 1e-15
# d3.options['N'] = 1
# d3.options['Area'] = 1

# d4 = y.add_diode('D4', 'gnd', 'n1')
# d4.options['Is'] = 1e-15
# d4.options['N'] = 1
# d4.options['Area'] = 1

# remember to always use the harmonics number to a power of 2
hb = HarmonicBalance2('HB1', 1e6, 10)

x = hb.run(y)

# print(x)
