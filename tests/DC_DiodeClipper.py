import logging
import numpy as np
import matplotlib.pyplot as plt

import setup
from Xyce import getXyceData
from yalrf import YalRF

y = YalRF("Circuit 3")

y.add_resistor('R1', 'n2', 'n3', 1e3)
y.add_resistor('R2', 'n1', 'n2', 3.3e3)
y.add_resistor('R3', 'n2', 'gnd', 3.3e3)
y.add_resistor('R4', 'n4', 'gnd', 5.6e3)
y.add_capacitor('C1', 'n2', 'n4', 0.47e-6)

y.add_vdc('Vcc', 'n1', 'gnd', 5)
vin = y.add_vdc('Vin', 'n3', 'gnd', 0)

d1 = y.add_diode('D1', 'n2', 'n1')
d2 = y.add_diode('D2', 'gnd', 'n2')

d1.options['Is'] = 4e-10
d1.options['Rs'] = 0 # not implemented yet
d1.options['N'] = 1.48
d1.options['Tt'] = 8e-7
d1.options['Cj0'] = 1.95e-11
d1.options['Vj'] = 0.4
d1.options['M'] = 0.38
d1.options['Eg'] = 1.36
d1.options['Xti'] = -8
d1.options['Kf'] = 0
d1.options['Af'] = 1
d1.options['Fc'] = 0.9
d1.options['Bv'] = 600
d1.options['Ibv'] = 1e-4
d2.options = d1.options.copy()

dc1 = y.add_dc_analysis('DC1')
#dc1.options['reltol'] = 1e-6

# get output from Xyce simulator
xyce = getXyceData('tests/data/circuit3.prn')

# sweep of Vin DC voltage
vsweep = np.arange(-10, 15, 1)
v2_yalrf = np.empty(len(vsweep))
v2_xyce = np.empty(len(vsweep))
i = 0
for v in vsweep:
    # run dc analysis with new dc voltage
    vin.dc = float(v)
    x = y.run('DC1')

    # get output data
    vn1 = y.get_voltage('DC1', 'n1')
    vn2 = y.get_voltage('DC1', 'n2')
    vn3 = y.get_voltage('DC1', 'n3')
    vn4 = y.get_voltage('DC1', 'n4')

    # compare outputs from Xyce and YalRF
    v2_yalrf[i] = vn2
    v2_xyce[i] = xyce[1][i][2]
    i = i + 1

plt.plot(vsweep, v2_yalrf)
plt.plot(vsweep, v2_xyce)
plt.title('Diode Clipper')
plt.grid()
plt.legend(['YalRF','Xyce'])
plt.xlabel('Vsweep')
plt.ylabel('Vclipped')
plt.show()