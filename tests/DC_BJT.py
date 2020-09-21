import numpy as np
import matplotlib.pyplot as plt

import setup
from Xyce import getXyceData
from yalrf import YalRF

y = YalRF("BJT Curve Tracer")

v1 = y.add_vdc('V1', 'nc', 'gnd', 1)
i1 = y.add_idc('I1', 'nb', 'gnd', 1e-6)

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

dc1 = y.add_dc_analysis('DC1')

vcesweep = np.linspace(0, 1, 101)
ibsweep = np.linspace(10e-6, 1010e-6, 11)

vbe_yalrf = []
vce_yalrf = []
ib_yalrf = []
ic_yalrf = []
for ib in ibsweep:
    vbe = np.zeros(len(vcesweep))
    ic = np.zeros(len(vcesweep))
    i = 0
    for vce in vcesweep:
        i1.dc = ib
        v1.dc = vce
        x = y.run('DC1')
        vbe[i] = y.get_voltage('DC1', 'nb')
        ic[i] = -x[2,0]
        i = i + 1
    vbe_yalrf.append(vbe)
    vce_yalrf.append(vcesweep)
    ib_yalrf.append(ibsweep)
    ic_yalrf.append(ic)

# get output from Xyce
xyce = np.array(getXyceData('tests/data/circuit4.prn')[1][:])

# number of sweep points in the base current
sweeplength = 11

# split the solution into arrays for each base current
vbe_xyce = np.split(xyce[:,1], sweeplength)
vce_xyce = np.split(xyce[:,2], sweeplength)
ib_xyce  = np.split(xyce[:,3], sweeplength)
ic_xyce  = np.split(-xyce[:,4], sweeplength)

plt.figure(figsize=(12,4))

plt.subplot(121)
for i in range(0, len(ibsweep)):
    plt.plot(vce_yalrf[i], ic_yalrf[i])
plt.title('YalRF BJT Curve Tracer')
plt.xlabel('Vce [V]')
plt.ylabel('Ic [mA]')
plt.grid()

plt.subplot(122)
for i in range(0, sweeplength):
    plt.plot(vce_xyce[i], ic_xyce[i])
plt.title('Xyce BJT Curve Tracer')
plt.xlabel('Vce [V]')
plt.ylabel('Ic [mA]')
plt.grid()

plt.show()