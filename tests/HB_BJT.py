import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance

y = YalRF('Diode Testbench')

i1 = y.add_iac('I1', 'n1', 'gnd', ac=10e-3, phase=+90)
i2 = y.add_idc('I2', 'n1', 'gnd', dc=10e-3)

i3 = y.add_idc('I3', 'n2', 'gnd', dc=20e-3)

r1 = y.add_resistor('R1', 'n1', 'nb', 1e3)
r2 = y.add_resistor('R2', 'nb', 'gnd', 1e3)

r3 = y.add_resistor('R3', 'n2', 'gnd', 100)
r4 = y.add_resistor('R4', 'n2', 'nc', 100)

q1 = y.add_bjt('Q1', 'nb', 'nc', 'gnd')

q1.options['Is'] = 1e-15
q1.options['Bf'] = 100
q1.options['Br'] = 1

hb = HarmonicBalance('HB1', 1e6, 5)

V = hb.run(y)
