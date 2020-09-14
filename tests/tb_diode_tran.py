import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF

y = YalRF('Diode Rectifier')

#y.add_vpulse('V1', 'n1', 'gnd', v1=0, v2=1, tstart=2e-9, tstop=6e-9, trise=1.5e-9, tfall=0.5e-9)
y.add_vsine('V1', 'n1', 'gnd', dc=0, ac=10, freq=60, phase=90)
y.add_resistor('R1', 'n2', 'gnd', 10e3)
y.add_capacitor('C1', 'n2', 'gnd', 10e-6)
#y.add_capacitor('C2', 'n1', 'n2', 1e-6)

d1 = y.add_diode('D1', 'n1', 'n2')
d1.options['Is'] = 4e-10
d1.options['N'] = 1.48
d1.options['Rs'] = 0.105
d1.options['Area'] = 1
d1.options['Cp'] = 1e-6
d1.options['Cj0'] = 1.95e-11
d1.options['Tt'] = 8e-7

tr1 = y.add_tran_analysis('TR1', tstop=100e-3, maxtstep=1e-3)
#tr1.options['max_iterations'] = 500
#tr1.options['reltol'] = 1e-6

x = y.run('TR1')

y.print_dc_voltages('TR1')
y.print_dc_currents('TR1')

t = y.get_time('TR1')
vn1 = y.get_voltage('TR1', 'n1')
vn2 = y.get_voltage('TR1', 'n2')
ir1 = y.get_itran('TR1', 'R1')
id1 = y.get_itran('TR1', 'D1')
iv1 = y.get_itran('TR1', 'V1')

#print(d1.Ic[:])
# print(d1.Id[:])
# print(vn1[:])
# print(vn2[:])
# print(iv1[:])

plt.figure(figsize=(12,5))

plt.subplot(121)
plt.plot(t, vn1)
plt.plot(t, vn2)
plt.title('Voltages')
plt.xlabel('Time [s]')
plt.ylabel('Voltage [V]')
plt.legend(['vn1', 'vn2', 'vn3'])
plt.grid()

plt.subplot(122)
plt.plot(t, ir1 * 1e3)
plt.plot(t, id1 * 1e3)
plt.title('Currents')
plt.xlabel('Time [s]')
plt.ylabel('Current [mA]')
plt.legend(['IR1', 'ID1'])
plt.grid()

plt.show()

