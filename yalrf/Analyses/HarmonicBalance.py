import numpy as np
from scipy import linalg

from yalrf.Netlist import Netlist
from yalrf.Analyses import AC, DC
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
                port = (dev.n1, dev.n2)
                if port not in nonlin_ports:
                    nonlin_ports.append(port)
            # elif isinstance(dev, BJT):
            #     nonlin_ports.append((dev.n1, dev.n2))
            #     nonlin_ports.append((dev.n1, dev.n3))

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
        
        T = 1. / self.freq          # period of the 1st harmonic
        dt = T / (4. * K)           # timestep in the transient waveform
        w1 = 2 * np.pi * self.freq  # angular frequency of the 1st harmonic

        # time array (oversampled 2x)
        time = np.linspace(0, T, 4 * K, endpoint=False)

        # array with the Kth harmonic frequencies
        freqs = self.freq * np.linspace(0, K, K+1)

        """ Transadmittance matrix (Y) """

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

        """ Interconnection current caused by voltage sources (Is) """

        Is = [np.zeros((M,1)) for i in range(K+1)]
        Is[0] = Yvs[0] @ Vs[0] # dc
        Is[1] = Yvs[1] @ Vs[1] # 1st harmonic

        """ Create initial voltage estimative for each port """

        # create port voltages using AC simulation results
        ac = AC('HB.AC', start=freqs[1], stop=freqs[K], numpts=K)
        ac.run(netlist)
        vac = ac.get_ac_solution()
        vdc = ac.get_dc_solution()

        # add dc component
        Vnode_td = [np.zeros((4 * K))]
        for n in range(L-1):
            Vnode_td.append(vdc[n] * np.ones((4 * K)))

        # add the ac components
        for n in range(L-1):
            for k in range(K):
                A = abs(vac[k,n])
                phi = (k + 1) * w1 * time + np.angle(vac[k,n])
                Vnode_td[n+1] = Vnode_td[n+1] + A * np.sin(phi)

        # calculate the initial port voltages (differential)
        Vport_td = []
        for port in nonlin_ports:
            Vport_td.append(Vnode_td[port[0]] - Vnode_td[port[1]])

        # obtain the spectrum of each port voltage (limit to K harmonics)
        Vport_fd = []
        for v in Vport_td:
            Vport_fd.append(scipy.fft.rfft(v)[0:K+1] / (4 * K))

        """ Time-domain g(t) waveforms """

        # create nonlinear subcircuit netlist
        nonlin_netlist = Netlist('Netlist of the nonlinear subciruit')
        for dev in nonlin_devs:
            if isinstance(dev, Diode):
                # get diode node names from full netlist
                n1_name = netlist.get_node_name(dev.n1)
                n2_name = netlist.get_node_name(dev.n2)

                # add the diode to the new netlist
                n1 = nonlin_netlist.add_node(n1_name)
                n2 = nonlin_netlist.add_node(n2_name)

                diode = Diode(dev.name, n1, n2)
                diode.options = dev.options
                nonlin_netlist.devices.append(diode)
            
            # TODO: generalize this for 3-port devices also

        # add a voltage source to each port of the nonlinear subcircuit
        vsources = []
        for port in nonlin_ports:
            n1_name = netlist.get_node_name(port[0])
            n2_name = netlist.get_node_name(port[1])
            
            vname = 'V_' + n1_name + '_' + n2_name
            v = nonlin_netlist.add_vdc(vname, n1_name, n2_name, 0)
            vsources.append(v)

        # get nonlinear subcircuits device list
        nldevs = nonlin_netlist.get_nonlinear_devices()

        # run a time-varying DC analysis to get the small-signal conductances over time
        gt = [[] for i in range(len(nldevs))]
        for i in range(len(time)):

            # set each voltage with its time-varying value
            for j in range(len(vsources)):
                vsources[j].dc = Vport_td[j][i]
            
            # run DC analysis
            dc = DC('HB.DC')
            dc.run(nonlin_netlist)

            # get the conductances from DC oppoint
            for k in range(len(nldevs)):
                gt[k].append(nldevs[k].gt)

        for dev in nldevs:
            if isinstance(dev, Diode):
                n1_name = nonlin_netlist.get_node_name(dev.n1)
                n2_name = nonlin_netlist.get_node_name(dev.n2)

                n1 = netlist.get_node_idx(n1_name)
                n2 = netlist.get_node_idx(n2_name)
                
                port = (n1, n2)
                print(nonlin_ports.index(port))

            # TODO: generalize this section for 3-port devices

        print(gt)
        # fourier transform of the port voltages

        # get the conductance waveforms g(t)

        # create the Jacobian and F(V0)

        # apply NR to obtain a new guess

        # fourier transform of the nonlinear currents

        # calculate F(V1) using the new voltage guess

        