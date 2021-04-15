import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import MultiToneHarmonicBalance

y = YalRF('Diode Testbench')

f = 10e6

i1 = y.add_iac('I1', 'n1', 'gnd', ac=10e-3, freq=f)
i2 = y.add_idc('I2', 'n1', 'gnd', dc=10e-3)

i3 = y.add_idc('I3', 'n2', 'gnd', dc=20e-3)

r1 = y.add_resistor('R1', 'n1', 'gnd', 1e3)
r2 = y.add_resistor('R2', 'n1', 'nb', 1e3)

r3 = y.add_resistor('R3', 'n2', 'gnd', 100)
r4 = y.add_resistor('R4', 'n2', 'nc', 1e3)

r5 = y.add_resistor('R5', 'ne', 'gnd', 1e3)

q1 = y.add_bjt('Q1', 'nb', 'nc', 'ne')

q1.options['Is'] = 1e-15
q1.options['Bf'] = 100
q1.options['Br'] = 1

# q1.options['Is'] = 1.4e-14
# q1.options['Nf'] = 1
# q1.options['Nr'] = 1
# q1.options['Ikf'] = 0.025
# q1.options['Ikr'] = 1e9
# q1.options['Vaf'] = 100
# q1.options['Var'] = 1e12
# q1.options['Ise'] = 3e-13
# q1.options['Ne'] = 1.5
# q1.options['Isc'] = 0
# q1.options['Nc'] = 2
# q1.options['Bf'] = 300
# q1.options['Br'] = 7.5

# q1.options['Cje'] = 4.5e-12
# q1.options['Vje'] = 0.75
# q1.options['Mje'] = 0.33
# q1.options['Cjc'] = 3.5e-12
# q1.options['Vjc'] = 0.75
# q1.options['Mjc'] = 0.33
# q1.options['Xcjc'] = 1
# q1.options['Cjs'] = 0
# q1.options['Vjs'] = 0.75
# q1.options['Mjs'] = 0
# q1.options['Fc'] = 0.5
# q1.options['Tf'] = 4e-10
# q1.options['Tr'] = 2.1e-8

hb = MultiToneHarmonicBalance('HB1', f, 10)

converged, freqs, Vf, time, Vt = hb.run(y)

hb.plot_v('nc')
plt.show()
