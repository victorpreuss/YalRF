import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import MultiToneHarmonicBalance

y = Netlist('Oscillator')

# circuit parameters
vcc = 5
c1 = 5e-9
c2 = 5e-9
c3 = 5e-9
r1 = 15e3
r2 = 15e3
r3 = 2.2e3
r4 = 5e3
re = 5

vin = 1
freq = 4e3

# VCC
y.add_idc('I1', 'nx', 'gnd', dc=vcc)
y.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

# vin oscprobe
Voscprobe = y.add_iac('I2', 'ny', 'gnd', ac=vin, phase=0)
y.add_gyrator('G2', 'ny', 'nz', 'gnd', 'gnd', 1)

# ideal harmonic filter of oscprobe
Zoscprobe = y.add_idealharmonicfilter('X1', 'nz', 'nc', freq)

# bias resistors
y.add_resistor('R4', 'nvcc', 'nc', r4)
y.add_resistor('RE', 'ne', 'gnd', re)

# RC feedback network
y.add_capacitor('C1', 'nc', 'n2', c1)
y.add_capacitor('C2', 'n2', 'nb', c2)
y.add_capacitor('C3', 'n1', 'gnd', c3)

y.add_resistor('R1', 'nc', 'n1', r1)
y.add_resistor('R2', 'n1', 'nb', r2)
y.add_resistor('R3', 'n2', 'gnd', r3)

# bjt
q1 = y.add_bjt('Q1', 'nb', 'nc', 'ne')

q1.options['Is'] = 10.3e-15
q1.options['Bf'] = 1090
q1.options['Br'] = 4
q1.options['Ne'] = 2
q1.options['Nc'] = 2
q1.options['Ise'] = 5.82e-12
q1.options['Vaf'] = 121
q1.options['Var'] = 24
q1.options['Ikf'] = 60.7e-3
q1.options['Ikr'] = 0.15

q1.options['Cje'] = 26.8e-12
q1.options['Vje'] = 1.1
q1.options['Mje'] = 0.5
q1.options['Cjc'] = 0 # 8.67e-12
q1.options['Vjc'] = 0.3
q1.options['Mjc'] = 0.3
q1.options['Xcjc'] = 1
q1.options['Cjs'] = 0
q1.options['Vjs'] = 0.75
q1.options['Mjs'] = 0
q1.options['Fc'] = 0.5
q1.options['Tf'] = 0 # 427e-12
q1.options['Tr'] = 0 # 50.3e-9

hb = MultiToneHarmonicBalance('HB1', freq, 10)
hb.options['maxiter'] = 100

Vprev = 0

def objFunc(x, info):

    global hb
    global Vprev

    # get new solution candidates
    fosc = x[0]
    Vosc = x[1]

    # update oscprobe values
    Voscprobe.ac = Vosc
    Voscprobe.freq = fosc
    Zoscprobe.freq = fosc

    # run harmonic balance

    hb.freq = fosc
    if info['itercnt'] > 0:
        converged, freqs, Vf, time, Vt = hb.run(y, Vprev)
    else:
        converged, freqs, Vf, time, Vt = hb.run(y)
    Vprev = hb.V

    # if HB failed to converge, return a bad convergence value to minimizer
    if not converged:
        return 1e2

    # get nodes of IdealHarmonicFilter
    n1 = hb.get_node_idx('nc')
    n2 = hb.get_node_idx('nz')

    # mag(Yosc) is the objective function to be minimized
    Voscx = Vf[n1,1]                             # voltage across oscprobe
    Iosc  = (Vf[n1,1] - Vf[n2,1]) * Zoscprobe.g  # current through oscprobe
    Yosc  = Iosc / Voscx

    info['itercnt'] += 1
    print('\nIter\tFreq [kHz]\tVosc [V]\tmag(Yosc)')
    print('{}\t{:.8f}\t{:.8f}\t{:.2e}\n'.format(info['itercnt'], fosc / 1e3, Vosc, abs(Yosc)))

    return abs(Yosc)

x0 = [ 3.8e3, 0.5]

result = optimize.fmin(func     = objFunc,
                       x0       = x0,
                       args     = ({'itercnt' : 0},),
                       xtol     = 1e-3,
                       ftol     = 1e-3,
                       maxfun   = 200,
                       disp     = True,
                       retall   = False,
                       full_output = True)

xopt = result[0]

# get solution and unnormalize it
fosc = xopt[0]
Vosc = xopt[1]

# update oscprobe values
Voscprobe.ac = Vosc
Voscprobe.freq = fosc
Zoscprobe.freq = fosc

# run harmonic balance
hb = MultiToneHarmonicBalance('HB1', fosc, 10)
converged, freqs, Vf, time, Vt = hb.run(y, Vprev)

print('Frequency of oscillation = {} Hz'.format(fosc))
print('Oscillation amplitude = {} V'.format(Vosc))

hb.plot_v('nc')
plt.show()

