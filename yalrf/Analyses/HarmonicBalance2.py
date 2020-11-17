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

class HarmonicBalance2:

    def __init__(self, name, freq, numharmonics=4):
        self.name = name
        self.freq = freq
        self.numharmonics = numharmonics

    def run(self, netlist):
        # get data from netlist
        devs = netlist.get_devices()
        lin_devs = netlist.get_linear_devices()
        nonlin_devs = netlist.get_nonlinear_devices()
        iidx = netlist.get_mna_extra_rows_dict()

        np.set_printoptions(precision=4, suppress=False, linewidth=150)

        """ Get nonlinear and voltage source ports """

        # get list of differential ports of the nonlinear subcircuit
        nonlin_ports = []
        for dev in nonlin_devs:
            if isinstance(dev, Diode):
                port = (dev.n1, dev.n2)
                if port not in nonlin_ports:
                    nonlin_ports.append(port)
            # elif isinstance(dev, BJT):
            #     nonlin_ports.append((dev.n1, dev.n3))
            #     nonlin_ports.append((dev.n2, dev.n3))

        # get list of differential ports of input voltage sources
        # also get voltage source values at DC and first harmonic
        Vs_dc = []
        Vs_ac = []
        vsource_ports = []
        for dev in lin_devs:
            if isinstance(dev, VoltageSource):
                port = (dev.n1, dev.n2)
                if port not in nonlin_ports:
                    vsource_ports.append(port)
                Vs_dc.append(dev.dc)
                Vs_ac.append(dev.ac * np.exp(1j * dev.phase))
        Vs = [np.array([Vs_dc]).T, np.array([Vs_ac]).T]

        print('Vs = {}'.format(Vs))

        """ Get basic information """

        K = self.numharmonics       # frequency of the 1st harmonic
        L = netlist.get_num_nodes() # total number of nodes in the netlist
        M = len(vsource_ports)      # number of voltage source ports
        N = len(nonlin_ports)       # number of nonlinear ports

        Kk = 2 * (K + 1)            # size of block matrices in the system
        W = Kk * N                  # total size of the matrix system 2 * (K+1) * N
        
        T = 1. / self.freq          # period of the 1st harmonic
        w1 = 2 * np.pi * self.freq  # angular frequency of the 1st harmonic

        S = 2 * (K + 1) - 1         # time array length for DFT/IDFT
        dt = T / S                  # timestep in the transient waveform

        # array with the Kth harmonic frequencies
        freqs = self.freq * np.linspace(0, K, K+1)

        """ Assembly of nonlinear subcircuit """

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
            Yvs.append(Ytrans[:N,N:]) # (NxM)
            Ynl.append(Ytrans[:N,:N]) # (NxN)

            # Yvs[-1] = Yvs[-1] + np.eye(len(Yvs[-1])) * 1e-12
            # Ynl[-1] = Ynl[-1] + np.eye(len(Ynl[-1])) * 1e-12

        # form the Y matrix based on Ynl
        Y = np.zeros((W, W))
        for i in range(N):
            for j in range(N):
                for k in range(K+1):
                    # get the Ymnk 2x2 matrix
                    Ymnk = np.zeros((2,2))
                    Ymnk[0,0] = +Ynl[k][i,j].real
                    Ymnk[0,1] = -Ynl[k][i,j].imag
                    Ymnk[1,0] = +Ynl[k][i,j].imag
                    Ymnk[1,1] = +Ynl[k][i,j].real

                    # put it in the Y matrix
                    I = Kk * i + 2 * k
                    J = Kk * j + 2 * k
                    Y[I:I+2,J:J+2] = Ymnk

        print('Y = {}'.format(Y))
                            
        """ Interconnection current created by voltage sources (Is) """

        Itempdc = Yvs[0] @ Vs[0] # dc
        Itempac = Yvs[1] @ Vs[1] # 1st harmonic

        Is = np.zeros((W, 1))
        for i in range(N):
            Is[Kk*i+0] = Itempdc[i].real
            Is[Kk*i+1] = Itempdc[i].imag
            Is[Kk*i+2] = Itempac[i].real / 2.
            Is[Kk*i+3] = Itempac[i].imag / 2.

        print('Is = {}'.format(Is))
        
        """ Initial voltage estimative for each port (V) """

        # create port voltages using AC simulation results
        ac = AC('HB.AC', start=freqs[1], stop=freqs[K], numpts=K)
        ac.run(netlist)
        vac = ac.get_ac_solution()
        vdc = ac.get_dc_solution()

        # add 'gnd' node
        vdc = np.insert(vdc, 0, 0, axis=0)
        vac = np.insert(vac, 0, 0, axis=1)

        V = np.zeros((W, 1))

        for i in range(N):
            port = nonlin_ports[i]

            vf = np.zeros(K+1, dtype=complex)
            vf[0] = (vdc[port[0]] - vdc[port[1]]) / 2.
            for k in range(1, K+1):
                vf[k] = (vac[k-1,port[0]] - vac[k-1,port[1]]) / 2.

            vt = np.zeros((S, 1))
            for s in range(S):
                vt[s] = vf[0].real
                for k in range(1, K+1):
                    vt[s] = vt[s] + 2 * (vf[k].real * np.cos(2. * np.pi * k * s / S) -
                                         vf[k].imag * np.sin(2. * np.pi * k * s / S))
                    vt[s] = min([vt[s], 0.7])
                    vt[s] = max([vt[s], -0.7])

            # plt.plot(vt)
            # plt.show()

            for i in range(N):
                for k in range(K+1):
                    vk = 0
                    for s in range(S):
                        vk = vk + vt[s] * (np.cos(2. * np.pi * k * s / S) -
                                           1j * np.sin(2. * np.pi * k * s / S))
                    vk = vk / S
                    V[Kk*i+2*k] = vk.real
                    V[Kk*i+2*k+1] = vk.imag

        # Forced good initial condition
        # V = np.array([[-3.0241e-01],
        #     [ 0.0000e+00],
        #     [ 6.2560e-04],
        #     [-7.3801e-01],
        #     [ 1.6367e-01],
        #     [-1.0291e-04],
        #     [ 8.3063e-04],
        #     [-5.9465e-02],
        #     [ 5.9560e-03],
        #     [-5.3677e-04],
        #     [ 9.1949e-04],
        #     [-2.1654e-02],
        #     [-8.2171e-03],
        #     [-1.3999e-03],
        #     [ 3.7401e-04],
        #     [-6.0768e-03],
        #     [-7.4935e-03],
        #     [-2.3699e-03],
        #     [-1.2872e-03],
        #     [ 5.4917e-04],
        #     [-4.1727e-03],
        #     [-2.5829e-03]])

        print('V = {}'.format(V))

        for a in range(15): # run 20 iterations of HB

            """ Time-domain v(t) """

            vt = np.zeros((N, S))
            for i in range(N):

                # assemble complex array of spectra for port 'i'
                vf = np.zeros(K+1, dtype=complex)
                for k in range(K+1):
                    vf[k] = V[Kk*i+2*k+0] + 1j * V[Kk*i+2*k+1]

                # compute inverse fourier transform of voltage waveform
                for s in range(S):
                    vt[i,s] = vf[0].real
                    for k in range(1, K+1):
                        vt[i,s] = vt[i,s] + 2 * (vf[k].real * np.cos(2. * np.pi * k * s / S) -
                                                 vf[k].imag * np.sin(2. * np.pi * k * s / S))
                
                # plt.plot(vt[i,:])
                # plt.show()

            """ Time-domain g(t) and i(t) waveforms """

            # run a time-varying DC analysis to get the small-signal conductances over time
            # TODO: wasting memory here (and also everywhere else)
            gt = np.zeros((N, N, S))
            it = np.zeros((N, N, S))
            for i in range(S):

                # set each voltage with its time-varying value
                for j in range(N):
                    vsources[j].dc = vt[j,i]
                    print('Vsource = {}'.format(vsources[j].dc))
                
                # run DC analysis
                dc = DC('HB.DC')
                dc.run(nonlin_netlist)

                # get the conductances from DC oppoint
                # TODO: generalize this section for 3-port devices
                for k in range(len(nldevs)):
                    dev = nldevs[k]
                    if isinstance(dev, Diode):
                        n1_name = nonlin_netlist.get_node_name(dev.n1)
                        n2_name = nonlin_netlist.get_node_name(dev.n2)

                        n1 = netlist.get_node_idx(n1_name)
                        n2 = netlist.get_node_idx(n2_name)

                        port = (n1, n2)
                        n = nonlin_ports.index(port)
                        
                        gt[n,n,i] = gt[n,n,i] + dev.oppoint['gd']
                        it[n,n,i] = it[n,n,i] + dev.oppoint['Id']

                        # print('gd = {}'.format(dev.oppoint['gd']))
                        # print('Id = {}'.format(dev.oppoint['Id']))

            print('vt = {}'.format(vt[0,:]))
            print('gt = {}'.format(gt[0,0,:]))
            print('it = {}'.format(it[0,0,:]))

            # compute the spectrum of each port conductance
            G = np.zeros((N, N, K+1), dtype=complex)
            for i in range(N):
                for k in range(K+1):
                    for s in range(S):
                        G[i,i,k] = G[i,i,k] + gt[i,i,s] * (np.cos(2. * np.pi * k * s / S) -
                                                           1j * np.sin(2. * np.pi * k * s / S))
                    G[i,i,k] = G[i,i,k] / S

            """ Nonlinear current (Ig) """

            # TODO: generalize this section for 3-port devices
            Ig = np.zeros((W, 1))
            for i in range(N):
                for k in range(K+1):
                    Igk = complex()
                    for s in range(S):
                        Igk = Igk + it[i,i,s] * (np.cos(2. * np.pi * k * s / S) -
                                                 1j * np.sin(2. * np.pi * k * s / S))
                    Igk = Igk / S
                    Ig[Kk*i+2*k] = Igk.real
                    Ig[Kk*i+2*k+1] = Igk.imag

            """ Creating the dIdV matrix from G(jw) """

            # creating matrices that form dI/dV
            Toe = np.zeros((W, W))
            Han = np.zeros((W, W))
            D = np.eye(W)

            # DC elements are 1/2 (from DFT definition used)
            for i in range(N):
                k = Kk * i
                D[k,k] = 0.5
                D[k+1,k+1] = 0.5

            # creating the Toeplitz matrix
            for i in range(N):
                for j in range(N):
                    Tmn = np.zeros((Kk, Kk))
                    for k in range(K+1):
                        for l in range(k, K+1):
                            
                            z = abs(k-l)
                            t = np.zeros((2,2))

                            t[0,0] = +G[i,j,z].real
                            t[0,1] = +G[i,j,z].imag
                            t[1,0] = -G[i,j,z].imag
                            t[1,1] = +G[i,j,z].real

                            Tmn[2*k:2*k+2,2*l:2*l+2] = t

                            t[0,1] = -t[0,1]
                            t[1,0] = -t[1,0]

                            Tmn[2*l:2*l+2,2*k:2*k+2] = t

                    # insert Tmn (Kk x Kk) submatrix into the proper position
                    I = Kk * i
                    J = Kk * j
                    Toe[I:I+Kk,J:J+Kk] = Tmn

            # creating the Hankel matrix
            for i in range(N):
                for j in range(N):
                    Hmn = np.zeros((Kk, Kk))
                    for k in range(K+1):
                        for l in range(K+1):
                            
                            z = k + l
                            h = np.zeros((2,2))

                            if z < K+1:
                                h[0,0] = +G[i,j,z].real
                                h[0,1] = +G[i,j,z].imag
                                h[1,0] = +G[i,j,z].imag
                                h[1,1] = -G[i,j,z].real
                            else:
                                z = 2 * K + 1 - z
                                h[0,0] = +G[i,j,z].real
                                h[0,1] = -G[i,j,z].imag
                                h[1,0] = -G[i,j,z].imag
                                h[1,1] = -G[i,j,z].real

                            Hmn[2*k:2*k+2,2*l:2*l+2] = h

                    # insert Hmn (Kk x Kk) submatrix into the proper position
                    I = Kk * i
                    J = Kk * j
                    Han[I:I+Kk,J:J+Kk] = Hmn

            print('V = {}'.format(V))

            dIdV = (Toe + Han) @ D
            J = Y + dIdV         # + Omega * dQdV
            F = Is + Y @ V + Ig  # + j * Omega * Q

            # print('Is = {}'.format(Is))
            # print('Y = {}'.format(Y))
            # print('Is + YV = {}'.format(Is + Y @ V))
            # # print('Ig = {}'.format(Ig))
            # print('F = {}'.format(F))
            # print('G = {}'.format(G))
            # print('Toe = {}'.format(Toe))
            # print('Han = {}'.format(Han))
            # print('D = {}'.format(D))
            # print('dIdV = {}'.format(dIdV))
            # print('J = {}'.format(J))

            V = V - linalg.inv(J) @ F
            print('V = {}'.format(V))

        # assemble complex array of spectra for port 'i'
        i = 0
        vf = np.zeros(K+1, dtype=complex)
        for k in range(K+1):
            vf[k] = V[Kk*i+2*k+0] + 1j * V[Kk*i+2*k+1]

        # compute inverse fourier transform of voltage waveform
        S = 4 * K
        vt = np.zeros(S)
        for s in range(S):
            vt[s] = vf[0].real
            for k in range(1, K+1):
                vt[s] = vt[s] + 2 * (vf[k].real * np.cos(2. * np.pi * k * s / S) -
                                     vf[k].imag * np.sin(2. * np.pi * k * s / S))
        

        plt.figure()
        plt.subplot(121)
        plt.plot(vt)
        plt.grid()
        plt.subplot(122)
        plt.stem(freqs, abs(vf), use_line_collection=True, markerfmt='^')
        for f,v in zip(freqs, vf):
            label = "({:.3f})".format(abs(v))
            plt.annotate(label, (f,abs(v)), textcoords="offset points", xytext=(0,10), ha='center')
        plt.grid()
        plt.show()
