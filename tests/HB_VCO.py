import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance, MultiToneHarmonicBalance

y = Netlist('Oscillator')

Cvaractor = 7.5e-12
osc_node = 'ne'

# VCC
y.add_idc('I1', 'nx1', 'gnd', dc=12)
y.add_gyrator('G1', 'nx1', 'nvcc', 'gnd', 'gnd', 1)

# vin oscprobe
Voscprobe = y.add_iac('I2', 'nx2', 'gnd', ac=0)
y.add_gyrator('G2', 'nx2', 'nz', 'gnd', 'gnd', 1)

# ideal harmonic filter of oscprobe
Zoscprobe = y.add_idealharmonicfilter('X1', 'nz', osc_node, 1)
Zoscprobe.g = 1e7

# passives
y.add_resistor('R1', 'nvcc', 'n1', 5.6e3)
y.add_resistor('R2', 'n1', 'gnd', 2.7e3)
y.add_resistor('R3', 'n4', 'gnd', 56)
y.add_resistor('Rload', 'nload', 'gnd', 50)

# in parallel to the resonating inductor
y.add_resistor('Rloss', 'n2', 'n3', 50e3)

y.add_inductor('L1', 'n1', 'nb', 100e-6)
y.add_inductor('L2', 'n2', 'n3', 490e-9)
y.add_inductor('L3', 'ne', 'n4', 100e-6)

y.add_capacitor('C1', 'nb', 'n2', 1000e-12)
y.add_capacitor('C2', 'n3', 'gnd', Cvaractor)
y.add_capacitor('C3', 'nb', 'ne', 4e-12)
y.add_capacitor('C4', 'ne', 'gnd', 8e-12)
y.add_capacitor('C5', 'ne', 'nload', 47e-12)

# bjts
q1 = y.add_bjt('Q1', 'nb', 'nvcc', 'ne')

q1.options['Is'] = 1.4e-14
q1.options['Nf'] = 1
q1.options['Nr'] = 1
q1.options['Ikf'] = 0.025
q1.options['Ikr'] = 1e9
q1.options['Vaf'] = 100
q1.options['Var'] = 1e12
q1.options['Ise'] = 3e-13
q1.options['Ne'] = 1.5
q1.options['Isc'] = 0
q1.options['Nc'] = 2
q1.options['Bf'] = 100
q1.options['Br'] = 7.5
q1.options['Cje'] = 4.5e-12
q1.options['Vje'] = 0.75
q1.options['Mje'] = 0.33
q1.options['Cjc'] = 3.5e-12
q1.options['Vjc'] = 0.75
q1.options['Mjc'] = 0.33
q1.options['Xcjc'] = 1
q1.options['Cjs'] = 0
q1.options['Vjs'] = 0.75
q1.options['Mjs'] = 0
q1.options['Fc'] = 0.5
q1.options['Tf'] = 4e-10
q1.options['Tr'] = 2.1e-8

hb = MultiToneHarmonicBalance('HB1', 1, 10)
hb.options['maxiter'] = 100
Vprev = 0

def objFunc(x, info):

    global osc_node
    global hb
    global Vprev

    # get new solution candidates
    fosc = x[0]
    Vosc = x[1]

    if fosc < 0:
        return 1e6

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

    # hb.print_v('ne')

    # get nodes of IdealHarmonicFilter
    n1 = hb.get_node_idx(osc_node)
    n2 = hb.get_node_idx('nz')

    # mag(Yosc) is the objective function to be minimized
    Voscx = Vf[n1,1]                             # voltage across oscprobe
    Iosc  = (Vf[n1,1] - Vf[n2,1]) * Zoscprobe.g  # current through oscprobe
    Yosc  = Iosc / Voscx

    info['itercnt'] += 1
    print('\nIter\tFreq [MHz]\tVosc [V]\tmag(Yosc)')
    print('{}\t{:.8f}\t{:.8f}\t{:.2e}\n'.format(info['itercnt'], fosc / 1e6, abs(Vosc), abs(Yosc)))

    return abs(Yosc)

x0 = [ 120e6, 0.1]
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
hb.freq = fosc

# run harmonic balance
converged, freqs, Vf, time, Vt = hb.run(y, Vprev)

print('Frequency of oscillation = {} Hz'.format(fosc))
print('Oscillation amplitude = {} V'.format(Vosc))

# hb.print_v('nload')
hb.print_v('ne')
# hb.print_v('nb')
# hb.plot_v('nload')
hb.plot_v('ne')
# hb.plot_v('nb')
plt.show()


