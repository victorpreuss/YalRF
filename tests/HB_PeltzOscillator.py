import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance, MultiToneHarmonicBalance


y = Netlist('Oscillator')

# circuit parameters
vcc = 10
r = 200e3
l = 0.5e-3
c = 10e-9
re = 50e3
freq = 0

# VCC
y.add_idc('I1', 'nx', 'gnd', dc=vcc)
y.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

# vin oscprobe
Voscprobe = y.add_iac('I2', 'ny', 'gnd', ac=0)
y.add_gyrator('G2', 'ny', 'nz', 'gnd', 'gnd', 1)

# ideal harmonic filter of oscprobe
Zoscprobe = y.add_idealharmonicfilter('X1', 'nz', 'nb', freq)

# tank circuit
y.add_resistor('R1', 'nvcc', 'nb', r)
y.add_inductor('L1', 'nvcc', 'nb', l)
y.add_capacitor('C1', 'nvcc', 'nb', c)

# emitter resistance
y.add_resistor('RE', 'ne', 'gnd', re)

# bjts
q1 = y.add_bjt('Q1', 'nb', 'nvcc', 'ne')
q2 = y.add_bjt('Q2', 'nvcc', 'nb', 'ne')

q1.options['Is'] = 1e-15
q1.options['Bf'] = 200
q1.options['Br'] = 1
q2.options = q1.options.copy()

hb = MultiToneHarmonicBalance('HB1', 1, 10)
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
        return 1e6

    # get nodes of IdealHarmonicFilter
    n1 = hb.get_node_idx('nb')
    n2 = hb.get_node_idx('nz')

    # mag(Yosc) is the objective function to be minimized
    Voscx = Vf[n1,1]                             # voltage across oscprobe
    Iosc  = (Vf[n1,1] - Vf[n2,1]) * Zoscprobe.g  # current through oscprobe
    Yosc  = Iosc / Voscx

    info['itercnt'] += 1
    print('\nIter\tFreq [kHz]\tVosc [V]\tmag(Yosc)')
    print('{}\t{:.8f}\t{:.8f}\t{:.2e}\n'.format(info['itercnt'], fosc / 1e3, Vosc, abs(Yosc)))

    return abs(Yosc)

b  = [(1.1e9, 1.3e9), (2, 3)]
x0 = [ 100e3, 0.5]
result = optimize.fmin(func     = objFunc,
                       x0       = x0,
                       args     = ({'itercnt' : 0},),
                       xtol     = 1e-5,
                       ftol     = 1e-5,
                       maxfun   = 100,
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

hb.plot_v('nb')
hb.plot_v('ne')
plt.show()


