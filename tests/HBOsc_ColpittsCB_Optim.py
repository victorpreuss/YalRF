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
vbias = 1
ibias = 1e-3
rl = 850
fosc = 50e6
l = 50e-9
ct = 1 / (l * (2 * np.pi * fosc) ** 2)

n = 0.2 # capactive divider ratio (feedback)
c1 = ct / (1 - n)
c2 = ct / n

numharmonics = 20

f0 = 1 / (2 * np.pi * np.sqrt(l * ct))
V0 = 2 * ibias * rl * (1 - n)

# declare harmonic balance solver
hb = MultiToneHarmonicBalance('HB1')
hb.options['maxiter'] = 100

# VCC
y.add_idc('I1', 'nx1', 'gnd', dc=vcc)
y.add_gyrator('G1', 'nx1', 'nvcc', 'gnd', 'gnd', 1)

# Vbias
y.add_idc('I2', 'ny1', 'gnd', dc=vbias)
y.add_gyrator('G2', 'ny1', 'nb', 'gnd', 'gnd', 1)

# Ibias
y.add_idc('I3', 'gnd', 'ne', dc=ibias)

# passives
y.add_resistor('Rl', 'nvcc', 'nc', rl)
y.add_inductor('L1', 'nvcc', 'nc', l)
C1 = y.add_capacitor('C1', 'nc', 'ne', c1)
C2 = y.add_capacitor('C2', 'ne', 'gnd', c2)

# bjts
q1 = y.add_bjt('Q1', 'nb', 'nc', 'ne')

q1.options['Is'] = 1e-15
q1.options['Bf'] = 100

# single iteration
converged, freqs, Vf, _, _ = hb.run_oscillator(y, f0, numharmonics, V0, 'nc')

hb.plot_v('nc')
hb.plot_v('ne')
plt.figure()
plt.plot(q1.Ic)
plt.show()

begin = time.time()

# loop varying 'n'
iter_cnt = 0
freqs = []
vtank = []
vtank2 = []
vtank3 = []
ipk = []
pwr = []
vtank_calc = []
eff = []
nvec = np.geomspace(0.05, 0.5, 20)
for n in nvec:
    # update netlist
    C1.C = ct / (1 - n)
    C2.C = ct / n
    
    # run oscillator analysis
    if iter_cnt == 0:
        hb.run_oscillator(y, f0, numharmonics, V0, 'nc')
    else:
        hb.run_oscillator(y, f0, numharmonics, V0, 'nc', useprev=True)

    # get results
    freqs.append(hb.freq)
    ipk.append(np.max(q1.Ic))
    vtank2.append(np.abs(hb.get_v('nc')[2]) / np.abs(hb.get_v('nc')[1]))
    vtank3.append(np.abs(hb.get_v('nc')[3]) / np.abs(hb.get_v('nc')[1]))
    vtank.append(np.abs(hb.get_v('nc')[1]))
    vtank_calc.append(2 * ibias * rl * (1 - n))
    eff.append(np.abs(hb.get_v('nc')[1])**2 / (2 * rl) / (vcc * ibias))

    iter_cnt += 1

end = time.time()

print('Shape of V: {}'.format(hb.V.shape))
print('Running time: {}'.format(end-begin))

plt.figure()
plt.plot(nvec, vtank, label='Simulated')
plt.plot(nvec, vtank_calc, label='Calculated')
plt.xlabel('$\eta$')
plt.ylabel('$V_{tank}$ [V]')
plt.grid()

plt.figure()
plt.plot(nvec, eff)
plt.xlabel('$\eta$')
plt.ylabel('Efficiency')
plt.grid()

plt.figure()
plt.plot(nvec, vtank2)
plt.xlabel('$\eta$')
plt.ylabel('$V_{tank2}$ [V]')
plt.grid()

plt.figure()
plt.plot(nvec, vtank3)
plt.xlabel('$\eta$')
plt.ylabel('$V_{tank3}$ [V]')
plt.grid()

plt.figure()
plt.plot(nvec, ipk)
plt.xlabel('$\eta$')
plt.ylabel('$I_{c,pk}$ [A]')
plt.grid()

plt.figure()
plt.plot(nvec, freqs)
plt.xlabel('$\eta$')
plt.ylabel('$f_{osc}$ [MHz]')
plt.grid()
plt.show()

# def objFunc(x, info):

#     global hb

#     # update netlist
#     C1.C = ct / (1 - x)
#     C2.C = ct / x

#     # initial condition for oscillator
#     f0 = 1 / (2 * np.pi * np.sqrt(l * ct))
#     V0 = 2 * ibias * rl * (1 - x)

#     # run oscillator analysis
#     hb.run_oscillator(y, f0, numharmonics, V0, 'nc')

#     info['itercnt'] += 1

#     return - (np.abs(hb.get_v('nc')[1])**2 / (2 * rl) / (vcc * ibias))

# begin = time.time()
# bounds = (0.05, 0.5)
# args = ({'itercnt' : 0},)
# xopt = optimize.minimize_scalar(fun      = objFunc,
#                                 bounds   = bounds,
#                                 method   = 'bounded',
#                                 args     = args,
#                                 tol      = 1e-3)

# end = time.time()

# print('Shape of V: {}'.format(hb.V.shape))
# print('Running time: {}'.format(end-begin))
# print(xopt)
# print(xopt.x)

# hb.print_v('nc')
# hb.plot_v('nc')
# plt.show()

