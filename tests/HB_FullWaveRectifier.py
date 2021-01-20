import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance, MultiToneHarmonicBalance

y = YalRF('Full-Wave Rectifier Testbench')

# i1 = y.add_iac('I1', 'vinp', 'gnd', ac=15e-3, phase=+90)
# i2 = y.add_idc('I2', 'vinp', 'gnd', dc=0e-3)

# r1 = y.add_resistor('R1', 'vinp', 'gnd', 100)
# r2 = y.add_resistor('R2', 'voutp', 'voutn', 100)

# d1 = y.add_diode('D1', 'vinp', 'voutp')
# d2 = y.add_diode('D2', 'gnd', 'voutp')
# d3 = y.add_diode('D3', 'voutn', 'vinp')
# d4 = y.add_diode('D4', 'voutn', 'gnd')

i1 = y.add_iac('I1', 'vinp', 'vinn', ac=15e-3, phase=+90, freq=1e6)
i2 = y.add_idc('I2', 'vinp', 'vinn', dc=5e-3)

r1 = y.add_resistor('R1', 'vinp', 'vinn', 100)
r2 = y.add_resistor('R2', 'voutp', 'gnd', 100)

d1 = y.add_diode('D1', 'vinp', 'voutp')
d2 = y.add_diode('D2', 'vinn', 'voutp')
d3 = y.add_diode('D3', 'gnd', 'vinp')
d4 = y.add_diode('D4', 'gnd', 'vinn')

d1.options['Is'] = 1e-15
d1.options['N'] = 1
d1.options['Area'] = 1

d2.options['Is'] = 1e-15
d2.options['N'] = 1
d2.options['Area'] = 1

d3.options['Is'] = 1e-15
d3.options['N'] = 1
d3.options['Area'] = 1

d4.options['Is'] = 1e-15
d4.options['N'] = 1
d4.options['Area'] = 1

# hb = HarmonicBalance('HB1', 1e6, 10)
hb = MultiToneHarmonicBalance('HB1', 1e6, 10)

converged, freqs, Vf, time, Vt = hb.run(y)

hb.plot_v('voutp')
plt.show()
