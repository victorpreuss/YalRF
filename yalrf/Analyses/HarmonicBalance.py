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

        """ Get basic lengths """

        K = self.numharmonics           # frequency of the 1st harmonic
        M = netlist.get_num_vsources()  # number of voltage sources
        N = netlist.get_num_nodes() - 1 # total number of nodes in the netlist without 'gnd'

        Kk = 2 * (K + 1)                # size of block matrices in the system
        W = Kk * N                      # total size of the matrix system 2 * (K+1) * N
        
        T = 1. / self.freq              # period of the 1st harmonic
        w1 = 2 * np.pi * self.freq      # angular frequency of the 1st harmonic

        S = 2 * (K + 1) - 1             # time array length for DFT/IDFT
        dt = T / S                      # timestep in the transient waveform

        # array with the Kth harmonic frequencies
        freqs = self.freq * np.linspace(0, K, K+1)

        """ Independent current sources (Is) """

        Is = np.zeros((W, 1))
        for dev in lin_devs:
            if isinstance(dev, CurrentSource):
                iac = dev.ac * np.exp(1j * dev.phase)
                n = dev.n1 - 1
                Is[Kk*n+0] = dev.dc
                Is[Kk*n+1] = 0.
                Is[Kk*n+2] = iac.real / 2.
                Is[Kk*n+3] = iac.imag / 2.

        print('Is = {}'.format(Is))

        """ Transadmittance matrix (Y) """

        Y = np.zeros((W, W))
        for k in range(K+1):

            # add ground node
            Yk = np.zeros((N+1, N+1), dtype=complex)

            # add currently supported devices stamps
            for dev in lin_devs:
                if isinstance(dev, (Resistor, Capacitor, Inductor)):
                    dev.add_ac_stamps(Yk, None, None, None, freqs[k])

            # remove ground node
            Yk = Yk[1:,1:]

            # populate the Y matrix using Yk
            for i in range(N):
                for j in range(N):
                    # get the Ymnk 2x2 matrix
                    Ymnk = np.zeros((2,2))
                    Ymnk[0,0] = +Yk[i,j].real
                    Ymnk[0,1] = -Yk[i,j].imag
                    Ymnk[1,0] = +Yk[i,j].imag
                    Ymnk[1,1] = +Yk[i,j].real

                    # put it in the Y matrix
                    I = Kk * i + 2 * k
                    J = Kk * j + 2 * k
                    Y[I:I+2,J:J+2] = Ymnk

        print('Y = {}'.format(Y))
        
        """ Initial voltage estimative for each port (V) """

        # create port voltages using AC simulation results
        ac = AC('HB.AC', start=freqs[1], stop=freqs[K], numpts=K)
        ac.run(netlist)
        vac = ac.get_ac_solution()
        vdc = ac.get_dc_solution()

        V = np.zeros((W, 1))
        for i in range(N):
            V[Kk*i] = vdc[i]
            for k in range(1, 2):#(1, K+1):
                vport = vac[k,i]
                V[Kk*i+2*k+0] = vport.real / 2.
                V[Kk*i+2*k+1] = vport.imag / 2.

        V = np.array([[-0.2],
                      [0],
                      [0],
                      [0.7],
                      [0.15],
                      [0],
                      [0],
                      [0.06],
                      [0],
                      [0],
                      [0],
                      [0]])

        print('V = {}'.format(V))

        for a in range(20): # run 10 iterations of HB

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

                # plt.plot(vt[0,:])
                # plt.show()

            """ Time-domain g(t) and i(t) waveforms """

            # run a time-varying DC analysis to get the small-signal conductances over time
            # TODO: wasting memory here (and also everywhere else)
            gt = np.zeros((N, N, S))
            it = np.zeros((N, N, S))

            dev = nonlin_devs[0]
            #dev.init()

            for s in range(S):
                for dev in nonlin_devs:
                    if isinstance(dev, Diode):
                        #dev.init()
                        dev.calc_oppoint(vt[:,s])

                        n1 = dev.n1-1
                        n2 = dev.n2-1

                        if n1 >= 0:
                            gt[n1,n1,s] = gt[n1,n1,s] + dev.oppoint['gd']
                            it[n1,n1,s] = it[n1,n1,s] + dev.oppoint['Id']

                        if n2 >= 0:
                            gt[n2,n2,s] = gt[n2,n2,s] + dev.oppoint['gd']
                            it[n2,n2,s] = it[n2,n2,s] - dev.oppoint['Id']
                        
            print('vt = {}'.format(vt))
            print('gt = {}'.format(gt))
            print('it = {}'.format(it))

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
                    Igk = 0 + 1j*0
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

            # DC elements are 1/2 (from DFT derivation)
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

            dIdV = (Toe + Han) @ D
            J = Y + dIdV         # + Omega * dQdV
            F = Is + Y @ V + Ig  # + j * Omega * Q

            for i in range(N):
                k = Kk * i
                J[k+1,k+1] = 1.

            # print('Is = {}'.format(Is))
            print('Y = {}'.format(Y))
            print('V = {}'.format(V))
            print('YV = {}'.format(Y @ V))
            print('Ig = {}'.format(Ig))
            print('F = {}'.format(F))
            # print('G = {}'.format(G))
            # print('Toe = {}'.format(Toe))
            # print('Han = {}'.format(Han))
            # print('D = {}'.format(D))
            # print('dIdV = {}'.format(dIdV))
            print('J = {}'.format(J))
            print('V = {}'.format(V))

            V = V - linalg.inv(J) @ F
            print('V = {}'.format(V))

        # assemble complex array of spectra for node 'i'
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

        vt_plot1 = vt.copy()
        vf_plot1 = vf.copy()
        
        plt.figure()
        plt.subplot(121)
        plt.plot(vt_plot1)
        plt.grid()
        plt.subplot(122)
        plt.stem(freqs, abs(vf_plot1), use_line_collection=True, markerfmt='^')
        for f, v in zip(freqs, vf_plot1):
            label = "({:.3f})".format(abs(v))
            plt.annotate(label, (f,abs(v)), textcoords="offset points", xytext=(0,10), ha='center')
        plt.grid()
        plt.show()


