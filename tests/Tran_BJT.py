import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF

y = YalRF('BJT AC Testbench')

y.add_vdc('V1', 'vcc', 'gnd', 5)
y.add_resistor('R1', 'vcc', 'nc', 200)
y.add_vsine('V2', 'nb', 'gnd', dc=0.65, ac=1e-3, freq=1e6, phase=0)

q1 = y.add_bjt('Q1', 'nb', 'nc', 'gnd')
q1.options['Is'] = 8.11e-14
q1.options['Nf'] = 1
q1.options['Nr'] = 1
q1.options['Ikf'] = 0.5
q1.options['Ikr'] = 0.225
q1.options['Vaf'] = 113
q1.options['Var'] = 24
q1.options['Ise'] = 1.06e-11
q1.options['Ne'] = 2
q1.options['Isc'] = 0
q1.options['Nc'] = 2
q1.options['Bf'] = 205
q1.options['Br'] = 4
q1.options['Cje'] = 0.
q1.options['Cjc'] = 0.
q1.options['Cjs'] = 0.

dc1 = y.add_dc_analysis('DC1')
dc1.options['max_iterations'] = 100
#dc1.options['vabstol'] = 1e-6
#dc1.options['reltol'] = 1e-3
xdc = y.run('DC1')

y.print_dc_voltages('DC1')
y.print_dc_currents('DC1')

tr1 = y.add_tran_analysis('TR1', tstop=10e-6, maxtstep=10e-9)
xtran = y.run('TR1', xdc)

time = y.get_time('TR1')
vout = y.get_voltage('TR1', 'nc')

plt.figure(figsize=(12,4))

plt.plot(time * 1e6, vout)
plt.title('BJT AC Testbench')
plt.xlabel('Time [us]')
plt.ylabel('Vout [V]')
plt.grid()

plt.show()


