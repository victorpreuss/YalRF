import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

import setup
from yalrf import YalRF, Netlist
from yalrf.Analyses import HarmonicBalance, MultiToneHarmonicBalance


y = Netlist('Oscillator')

# circuit parameters
vcc = 5
c1 = 3.3e-12
c2 = 3.9e-12
r1 = 3e3
r2 = 6.8e3
re = 1.5e3
l1 = 6.944e-9

vin = 0
freq = 1.2e9

# VCC
y.add_idc('I1', 'nx', 'gnd', dc=vcc)
y.add_gyrator('G1', 'nx', 'nvcc', 'gnd', 'gnd', 1)

# vin oscprobe
Voscprobe = y.add_iac('I2', 'ny', 'gnd', ac=vin, phase=0)
y.add_gyrator('G2', 'ny', 'nz', 'gnd', 'gnd', 1)

# ideal harmonic filter of oscprobe
Zoscprobe = y.add_idealharmonicfilter('X1', 'nz', 'nind', freq)

# bias resistors
y.add_resistor('R1', 'nvcc', 'nb', r1)
y.add_resistor('R2', 'nb', 'gnd', r2)

# emitter resistance
y.add_resistor('RE', 'ne', 'gnd', re)

# capacitor feedback network
y.add_capacitor('C1', 'nb', 'ne', c1)
y.add_capacitor('C2', 'ne', 'gnd', c2)

# resonating inductor
y.add_inductor('L1', 'nind', 'gnd', l1)

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

hb = MultiToneHarmonicBalance('HB1', 1, 100)
hb.options['maxiter'] = 100

global result
global Vprev

def objFunc(x, info):

    global hb
    global Vprev
    global result

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
    n1 = hb.get_node_idx('nind')
    n2 = hb.get_node_idx('nz')

    # mag(Yosc) is the objective function to be minimized
    Voscx = Vf[n1,1]                             # voltage across oscprobe
    Iosc  = (Vf[n1,1] - Vf[n2,1]) * Zoscprobe.g  # current through oscprobe
    Yosc  = Iosc / Voscx

    info['itercnt'] += 1
    print('\nIter\tFreq [GHz]\tVosc [V]\tmag(Yosc)')
    print('{}\t{:.8f}\t{:.8f}\t{:.2e}\n'.format(info['itercnt'], fosc / 1e9, Vosc, abs(Yosc)))

    if abs(Yosc) < 1e-12:
        result = x
        print("Optimization terminated succesfully! |Y| < 1e-12")
        raise SmallEnoughException()

    return abs(Yosc)

class SmallEnoughException(Exception):
    pass

# result = optimize.fmin_bfgs(f        = objFunc,
#                             x0       = x0,
#                             args     = ({'itercnt' : 0},),
#                             gtol     = 1e-5,
#                             epsilon  = 1e-7,
#                             maxiter  = 10,
#                             disp     = True,
#                             retall   = False,
#                             full_output = True)[0]

# b  = [(1.1e9, 1.3e9), (2, 3)]
# result = optimize.minimize(fun=objFunc, x0=x0, method='Powell', bounds=b, args=({'itercnt' : 0},))[0]

x0 = [ 1.2e9, 1]
try:
    xopt = optimize.fmin(func     = objFunc,
                         x0       = x0,
                         args     = ({'itercnt' : 0},),
                         xtol     = 1e-3,
                         ftol     = 1e-3,
                         maxfun   = 200,
                         disp     = True,
                         retall   = False,
                         full_output = True)[0]
except SmallEnoughException:
    xopt = result


# get solution
fosc = xopt[0]
Vosc = xopt[1]

# update oscprobe values
Voscprobe.ac = Vosc
Voscprobe.freq = fosc
Zoscprobe.freq = fosc

# run harmonic balance
hb = MultiToneHarmonicBalance('HB1', fosc, 100)
converged, freqs, Vf, time, Vt = hb.run(y, Vprev)

print('Frequency of oscillation = {} Hz'.format(fosc))
print('Oscillation amplitude = {} V'.format(Vosc))

# hb.plot_v('nind')
# hb.plot_v('nl')
plt.show()


