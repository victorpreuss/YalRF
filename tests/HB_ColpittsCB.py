import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance, MultiToneHarmonicBalance

y = Netlist('Oscillator')

# circuit parameters
vcc = 3
vee = -3
rb = 100
re = 2.2e3
l = 5e-6
c1 = 200e-12
c2 = 200e-12
freq = 0

# VCC
y.add_idc('I1', 'nx1', 'gnd', dc=vcc)
y.add_gyrator('G1', 'nx1', 'nvcc', 'gnd', 'gnd', 1)

# VEE
y.add_idc('I2', 'nx2', 'gnd', dc=vee)
y.add_gyrator('G2', 'nx2', 'nvee', 'gnd', 'gnd', 1)

# vin oscprobe
Voscprobe = y.add_iac('I3', 'nx3', 'gnd', ac=0)
y.add_gyrator('G3', 'nx3', 'nz', 'gnd', 'gnd', 1)

# ideal harmonic filter of oscprobe
Zoscprobe = y.add_idealharmonicfilter('X1', 'nz', 'nc', freq)
Zoscprobe.g = 1e6

# passives
y.add_resistor('Rb', 'nb', 'gnd', rb)
y.add_resistor('Re', 'ne', 'nvee', re)
y.add_resistor('Rx', 'nvcc', 'nc', 2e3)
y.add_inductor('L1', 'nvcc', 'nc', l)
y.add_capacitor('C1', 'nvcc', 'ne', c1)
y.add_capacitor('C2', 'ne', 'nc', c2)

# bjts
q1 = y.add_bjt('Q1', 'nb', 'nc', 'ne')

q1.options['Is'] = 1e-16
q1.options['Bf'] = 100
q1.options['Br'] = 1

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

x0 = [ 5e6, .5]
result = optimize.fmin(func     = objFunc,
                       x0       = x0,
                       args     = ({'itercnt' : 0},),
                       xtol     = 1e-5,
                       ftol     = 1e-5,
                       maxfun   = 150,
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

hb.print_v('nc')
hb.plot_v('nc')
plt.show()


