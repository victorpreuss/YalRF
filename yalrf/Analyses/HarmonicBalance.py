import numpy as np
from scipy import linalg

from yalrf.Devices import *
from yalrf.Utils import hb_logger as logger

import logging
logger.setLevel(logging.INFO)

import scipy.fft
import matplotlib.pyplot as plt

class HarmonicBalance:

    def __init__(self, name, freq, numharmonics=3):
        self.name = name
        self.freq = freq
        self.numharmonics = numharmonics

    def run(self, netlist):
        # get data from netlist
        devs = netlist.get_devices()
        lin_devs = netlist.get_linear_devices()
        nonlin_devs = netlist.get_nonlinear_devices()
        iidx = netlist.get_mna_extra_rows_dict()

        # get list of differential ports of the nonlinear subcircuit
        nonlin_ports = []
        for dev in nonlin_devs:
            if isinstance(dev, Diode):
                nonlin_ports.append((dev.n1, dev.n2)) # anode - cathode
            elif isinstance(dev, BJT):
                nonlin_ports.append((dev.n1, dev.n3)) # base - emitter
                nonlin_ports.append((dev.n2, dev.n3)) # collector - emitter

        # get list of differential ports of input voltage sources
        # also get voltage source values at DC and first harmonic
        Vs_dc = []
        Vs_ac = []
        vsource_ports = []
        for dev in lin_devs:
            if isinstance(dev, VoltageSource):
                vsource_ports.append((dev.n1, dev.n2))
                Vs_dc.append(dev.dc)
                Vs_ac.append(dev.ac)
        Vs = [np.array([Vs_dc]).T, np.array([Vs_ac]).T]

        # get the appropriate lengths
        K = self.numharmonics
        L = netlist.get_num_nodes()
        M = len(vsource_ports)
        N = len(nonlin_ports)
        
        # period of the 1st harmonic
        T = 1. / self.freq
        dt = T / (2. * K)

        # time array
        time = np.linspace(0, T, 2 * K, endpoint=False)

        # array with the Kth harmonic frequencies
        freqs = self.freq * np.linspace(0, K, K+1)

        # create MNA matrices for linear subcircuit
        A = np.zeros((L, L), dtype=complex)
        C = np.zeros((N+M, L), dtype=complex)

        Yvs = [] # relates voltage source to interconnection current
        Ynl = [] # relates interconnection voltage to interconnection current
        for freq in freqs:
            # refresh matrices
            A[:,:] = 0.
            C[:,:] = 0.

            # add currently supported devices stamps
            for dev in lin_devs:
                if isinstance(dev, (Resistor, Capacitor, Inductor)):
                    dev.add_ac_stamps(A, None, None, None, freq)

            k = 0
            for node in nonlin_ports + vsource_ports:
                # get nodes
                n1 = node[0]
                n2 = node[1]

                # add to connexion matrix 
                C[k,n1] = +1.
                C[k,n2] = -1.

                # add 100 ohms resistor in parallel with nonlinear or vsource devices
                g = 0.01 
                A[n1][n1] = A[n1][n1] + g
                A[n2][n2] = A[n2][n2] + g
                A[n1][n2] = A[n1][n2] - g
                A[n2][n1] = A[n2][n1] - g
                
                k = k + 1

            # remove ground nodes
            An = A[1:,1:]
            Cn = C[:,1:]

            # calculate transimpedance and transdmittance matrices
            Ztrans = Cn @ linalg.pinv(An) @ Cn.transpose()
            Ytrans = linalg.pinv(Ztrans) - np.eye(N+M) * 0.01

            # split the transadmittance matrix between nonlinear and voltage source ports
            Yvs.append(Ytrans[:N,N:])
            Ynl.append(Ytrans[:N,:N])

        # calculate interconnection current caused by voltage sources
        Is = [np.zeros((M,1)) for i in range(0, K+1)]
        Is[0] = Yvs[0] @ Vs[0] # dc
        Is[1] = Yvs[1] @ Vs[1] # 1st harmonic

        # create initial voltage estimative for the ports
        vt_nodes = []
        vt_nodes.append(np.zeros((1, len(time))))
        for node in range(1, L):
            omega = 2 * np.pi * self.freq
            rand_amp = np.random.random_sample()
            rand_ph = 2. * np.pi * np.random.random_sample()
            rand_vt = rand_amp * np.sin(omega * time + rand_ph)
            vt_nodes.append(rand_vt)

        print(vt_nodes)

        vt_ports = []
        for node in nonlin_ports:
            n1 = node[0]
            n2 = node[1]
            vt_ports.append(vt_nodes[n1] - vt_nodes[n2])

        
        a = scipy.fft.fft(vt_ports[0])
        a = scipy.fft.fftshift(a)

        #a = np.fft.fft(vt_ports[0])
        #a = np.fft.fftshift(a)
        
        #idx = np.fft.fftfreq(len(a), dt)
        #idx = np.fft.fftshift(idx)

        plt.figure()
        plt.plot(abs(a))
        plt.grid()
        plt.show()

        # fourier transform of the port voltages

        # get the conductance waveforms g(t)

        # create the Jacobian and F(V0)

        # apply NR to obtain a new guess

        # fourier transform of the nonlinear currents

        # calculate F(V1) using the new voltage guess

        