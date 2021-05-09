import time
import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist

y = YalRF('Colpitts CC Oscillator')

# circuit parameters
vcc = 5
c1 = 3.3e-12
c2 = 3.9e-12
r1 = 3e3
r2 = 6.8e3
re = 1.5e3
l1 = 6.944e-9

# VCC
v1 = y.add_vdc('V1', 'nvcc', 'gnd', vcc)

# AC supply
i1 = y.add_iac('I1', 'nind', 'gnd', 1)

# bias resistors
y.add_resistor('R1', 'nvcc', 'nb', r1)
y.add_resistor('R2', 'nb', 'gnd', r2)

# emitter resistance
y.add_resistor('RE', 'ne', 'gnd', re)

# capacitor feedback network
y.add_capacitor('C1', 'nb', 'ne', c1)
y.add_capacitor('C2', 'ne', 'gnd', c2)

# dc feed and dc block (TODO: improve)
y.add_inductor('Lfeed', 'nvcc', 'nc', 1e-3)
y.add_capacitor('Cblk1', 'nb', 'nind', 1e-6)

# load connection
y.add_capacitor('Cblk2', 'nc', 'nl', 1e-6)
y.add_resistor('RL', 'nl', 'gnd', 50)

# bjt
q1 = y.add_bjt('Q1', 'nb', 'nc', 'ne')

q1.options['Is'] = 1e-15
q1.options['Bf'] = 100
q1.options['Br'] = 5
q1.options['Vaf'] = 60
q1.options['Var'] = 20

# run DC sim
dc1 = y.add_dc_analysis('DC1')
xdc = y.run('DC1')

y.print_dc_voltages('DC1')

# run AC sim
ac1 = y.add_ac_analysis('AC1', start=1e9, stop=1.8e9, numpts=500, sweeptype='linear')
xac = y.run('AC1', xdc)

freqs = y.get_freqs('AC1') / 1e9
vin = y.get_voltage('AC1', 'nind')
iin = i1.ac

Rin = np.real(vin / iin)
Xin = np.imag(vin / iin)

# point of interest ~1.4GHz
freq = freqs[265]
rin = Rin[265]
xin = Xin[265]
print('{:.2f} {:.2f} {:.2f}'.format(freq, rin, xin))

plt.figure(figsize=(12,4))

plt.subplot(121)
plt.plot(freqs, Rin)
plt.plot(freq, rin, marker='o')
plt.text(freq, 0.97*rin, '{:.2f}'.format(rin), ha='right', fontsize=11)
# plt.title('Input Resistance')
plt.xlabel('Frequency [GHz]')
plt.ylabel('Rin [Ohms]')
plt.grid()

plt.subplot(122)
plt.plot(freqs, Xin)
plt.plot(freq, xin, marker='o')
plt.text(freq, 0.97*xin, '{:.2f}'.format(xin), ha='right', fontsize=11)
# plt.title('Input Reactance')
plt.xlabel('Frequency [GHz]')
plt.ylabel('Xin [Ohms]')
plt.grid()
plt.show()


