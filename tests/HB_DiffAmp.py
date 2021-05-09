import numpy as np
import matplotlib.pyplot as plt

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import AC, HarmonicBalance, MultiToneHarmonicBalance

y = YalRF('Differential Amplifier')

# circuit parameters
ibias = 50e-3
vbias = 1.2
vcc = 5
vin = 60e-3
rload = 50

# VCC
i1 = y.add_idc('I1', 'nx', 'gnd', dc=vcc)
g1 = y.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

# vin1
i2 = y.add_iac('I2', 'ny', 'gnd', ac=vin, phase=0, freq=10e6)
i3 = y.add_idc('I3', 'ny', 'gnd', dc=vbias)
g2 = y.add_gyrator('G2', 'ny', 'nb1', 'gnd', 'gnd', 1)

# vin2
i4 = y.add_iac('I4', 'nz', 'gnd', ac=vin, phase=+180, freq=10e6)
i5 = y.add_idc('I5', 'nz', 'gnd', dc=vbias)
g3 = y.add_gyrator('G3', 'nz', 'nb2', 'gnd', 'gnd', 1)

# collector loads
r1 = y.add_resistor('R1', 'nvcc', 'nc1', rload)
r2 = y.add_resistor('R2', 'nvcc', 'nc2', rload)

# differential pair
q1 = y.add_bjt('Q1', 'nb1', 'nc1', 'ne')
q2 = y.add_bjt('Q2', 'nb2', 'nc2', 'ne')

q3 = y.add_bjt('Q3', 'nbx', 'ne', 'gnd')
q4 = y.add_bjt('Q4', 'nbx', 'nbx', 'gnd')

i6 = y.add_idc('I6', 'nbx', 'gnd', dc=ibias)

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
q1.options['Cje'] = 2.95e-11
q1.options['Cjc'] = 1.52e-11
q1.options['Cjs'] = 0.
# q1.options['Is'] = 1e-15
# q1.options['Bf'] = 100
# q1.options['Vaf'] = 60

q2.options = q1.options
q3.options = q1.options
q4.options = q1.options

# run harmonic balance
# hb = HarmonicBalance('HB1', 10e6, 10)
hb = MultiToneHarmonicBalance('HB1', 10e6, 20)
# hb.options['maxiter'] = 100
# converged, freqs, Vf, time, Vt = hb.run(y)

# hb.print_v('nc1')
# hb.plot_v('nc1')
# hb.plot_v('nc2')
# plt.show()

vi = []
vout = []
V0 = None
for vin in np.arange(100e-6, 50.1e-3, 2e-3):

    # update input voltage
    i2.ac = vin
    i4.ac = vin

    # run harmonic balance
    converged, freqs, Vf, time, Vt = hb.run(y, V0)
    V0 = hb.V
    
    # get input and output information
    vi.append(hb.get_v('nb1')[1] - hb.get_v('nb2')[1])
    vout.append(hb.get_v('nc2')[1])
    
vi = np.array(vi, dtype=complex)
vout = np.array(vout, dtype=complex)
gain = 20 * np.log10(np.abs(vout / vi))

ac1 = y.add_ac_analysis('AC1', start=1e6, stop=10e9, numpts=1000, sweeptype='linear')
xac = y.run('AC1')

freqs = y.get_freqs('AC1')
vi_ac = y.get_voltage('AC1', 'nb1') - y.get_voltage('AC1', 'nb2')
vout_ac = y.get_voltage('AC1', 'nc2')
gain_ac = 20 * np.log10(np.abs(vout_ac / vi_ac))

plt.figure(figsize=(14,5))

vv = np.abs(vi) * 1e3
plt.subplot(121)
plt.plot(vv, gain)
plt.xlabel('Input Voltage [mV]')
plt.ylabel('Voltage Gain [dB]')
plt.title('HB Gain Compression')
plt.grid()

idx = np.argmin(np.abs(gain-gain[0]+1))
plt.plot(vv[0], gain[0], color='red', marker='o')
label = '{:.1f} dB @ \n{:.0f} uV'.format(gain[0], vv[0] * 1e3)
plt.annotate(label, (vv[0] , gain[0]), textcoords="offset points", xytext=(5,-25), ha='left', va='bottom', rotation=0)
plt.plot(vv[idx], gain[idx], color='red', marker='o')
label = '{:.1f} dB @ \n{:.0f} mV'.format(gain[idx], vv[idx])
plt.annotate(label, (vv[idx] , gain[idx]), textcoords="offset points", xytext=(0,5), ha='left', va='bottom', rotation=0)

plt.subplot(122)
plt.semilogx(freqs, gain_ac, color='red')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Voltage Gain [dB]')
plt.title('AC Small-Signal Gain')
plt.grid()

idx = np.argmin(np.abs(gain_ac-gain_ac[0]+3))
plt.plot(freqs[0], gain_ac[0], color='blue', marker='o')
label = '{:.1f} dB @ \n{:.0f} MHz'.format(gain_ac[0], freqs[0] / 1e6)
plt.annotate(label, (freqs[0] , gain_ac[0]), textcoords="offset points", xytext=(5,-25), ha='left', va='bottom', rotation=0)
plt.plot(freqs[idx], gain_ac[idx], color='blue', marker='o')
label = '{:.1f} dB @ \n{:.0f} MHz'.format(gain_ac[idx], freqs[idx] / 1e6)
plt.annotate(label, (freqs[idx] , gain_ac[idx]), textcoords="offset points", xytext=(5,0), ha='left', va='bottom', rotation=0)

plt.show()




