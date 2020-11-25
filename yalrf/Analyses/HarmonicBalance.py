import numpy as np
from scipy import linalg

from yalrf.Netlist import Netlist
from yalrf.Analyses import AC, DC
from yalrf.Devices import *
from yalrf.Utils import hb_logger as logger

import logging
logger.setLevel(logging.INFO)

import matplotlib.pyplot as plt

options = dict()
options['reltol'] = 1e-3
options['abstol'] = 1e-6
options['maxiter'] = 100

class HarmonicBalance:

    def __init__(self, name, freq, numharmonics=4):
        self.name = name
        self.freq = freq
        self.numharmonics = numharmonics

        self.options = options.copy()

    def run(self, netlist):

        self.netlist = netlist.copy()

        # get device data from netlist
        self.lin_devs = self.netlist.get_linear_devices()
        self.nonlin_devs = self.netlist.get_nonlinear_devices()

        # print option to make large matrices readable
        np.set_printoptions(precision=4, suppress=False, linewidth=150)

        """ Get basic sizes """

        K = self.numharmonics           # frequency of the 1st harmonic
        N = netlist.get_num_nodes() - 1 # total number of nodes in the netlist without 'gnd'

        Kk = 2 * (K + 1)                # size of block matrices in the system
        W = Kk * N                      # total size of the matrix system 2 * (K+1) * N
        
        T = 1. / self.freq              # period of the 1st harmonic
        S = 2 * K + 1                   # time array length for DFT/IDFT

        # array with the Kth harmonic frequencies
        freqs = self.freq * np.linspace(0, K, K+1)

        self.W = W
        self.Kk = Kk
        self.N = N
        self.K = K
        self.S = S
        self.freqs = freqs

        """ Independent current sources (Is) """

        Is = self.calc_Is()

        print('Is = {}'.format(Is))

        """ Transadmittance matrix (Y) """

        Y = self.calc_Y()

        print('Y = {}'.format(Y))
        
        """ Initial voltage estimation for each node (V, vt) """

        V, vt = self.calc_V0()

        print('V = {}'.format(V))

        """ Get simulation options """

        abstol = self.options['abstol']
        reltol = self.options['reltol']
        maxiter = self.options['maxiter']

        """ Start Harmonic Balance Solver """

        gt = np.zeros((N+1, N+1, S))
        G = np.zeros((self.N, self.N, self.K+1), dtype=complex)

        it = np.zeros((N+1, S))
        Inl = np.zeros((W, 1))

        itercnt = 0
        while True:

            """ Time-domain v(t) """

            vtold = vt.copy()

            self.ifft_V(V, vt)

            # print(vt)
            # plt.plot(vt[0,:])
            # plt.show()

            """ Time-domain g(t) and i(t) waveforms """

            # run a time-varying oppoint analysis on the nonlinear devices
            gt[:,:,:] = 0
            it[:,:] = 0
            for s in range(self.S):
                for dev in self.nonlin_devs:
                    dev.add_hb_stamps(vt, it, gt, s)

            # print('vt = {}'.format(vt))
            # print('gt = {}'.format(gt))
            # print('it = {}'.format(it))
            
            """ Nonlinear conductance G(jw) """

            self.fft_G(gt[1:,1:,:], G)

            """ Nonlinear current Inl(jw) """

            self.fft_Inl(it[1:,:], Inl)

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
            J = Y + dIdV                    # + Omega * dQdV

            for i in range(N):
                k = Kk * i
                J[k+1,k+1] = 1.

            Il = Y @ V - Is

            converged = True
            
            F = Il + Inl                    # + j * Omega * Q
            for i in range(W):
                if abs(F[i]) > abstol:
                    converged = False

            Frel = 2 * abs(F / (Il - Inl + 1e-12))
            for i in range(W):
                if Frel[i] > reltol:
                    converged = False

            itercnt += 1
            if converged or itercnt >= maxiter:
                break

            V = V - linalg.inv(J) @ F

            # print('Is = {}'.format(Is))
            # print('Y = {}'.format(Y))
            # print('YV = {}'.format(Y @ V))
            # print('Ig = {}'.format(Ig))
            # print('F = {}'.format(F))
            # print('G = {}'.format(G))
            # print('Toe = {}'.format(Toe))
            # print('Han = {}'.format(Han))
            # print('D = {}'.format(D))
            # print('dIdV = {}'.format(dIdV))
            # print('J = {}'.format(J))
            # print('V = {}'.format(V))

        S = 8 * K   # increase number of time samples
        Vt = np.zeros((N, S))
        Vf = np.zeros((N, K+1), dtype=complex)
        for i in range(N):

            # assemble complex array of spectra for node 'i'
            for k in range(K+1):
                Vf[i,k] = V[Kk*i+2*k+0] + 1j * V[Kk*i+2*k+1]

            # compute inverse fourier transform of voltage waveform
            for s in range(S):
                Vt[i,s] = Vf[i,0].real
                for k in range(1, K+1):
                    Vt[i,s] = Vt[i,s] + 2 * (Vf[i,k].real * np.cos(2. * np.pi * k * s / S) -
                                             Vf[i,k].imag * np.sin(2. * np.pi * k * s / S))

        return freqs, Vf, Vt

    def calc_Is(self):
        Kk = self.Kk
        Is = np.zeros((self.W, 1))
        for dev in self.lin_devs:
            if isinstance(dev, CurrentSource):
                iac = dev.ac * np.exp(1j * dev.phase)

                n1 = dev.n1 - 1
                if n1 >= 0:
                    Is[Kk*n1+0] += dev.dc
                    Is[Kk*n1+1] += 0.
                    Is[Kk*n1+2] += iac.real / 2.
                    Is[Kk*n1+3] += iac.imag / 2.

                n2 = dev.n2 - 1
                if n2 >= 0:
                    Is[Kk*n2+0] -= dev.dc
                    Is[Kk*n2+1] -= 0.
                    Is[Kk*n2+2] -= iac.real / 2.
                    Is[Kk*n2+3] -= iac.imag / 2.

        return Is

    def calc_Y(self):

        Y = np.zeros((self.W, self.W))
        for k in range(self.K+1):

            # add ground node
            Yk = np.zeros((self.N+1, self.N+1), dtype=complex)

            # add currently supported devices stamps
            for dev in self.lin_devs:
                if isinstance(dev, (Resistor, Capacitor, Inductor, Gyrator)):
                    dev.add_ac_stamps(Yk, None, None, None, self.freqs[k])

            # remove ground node
            Yk = Yk[1:,1:]

            # populate the Y matrix using Yk
            for i in range(self.N):
                for j in range(self.N):

                    # get the Ymnk 2x2 matrix
                    Ymnk = np.zeros((2,2))
                    Ymnk[0,0] = +Yk[i,j].real
                    Ymnk[0,1] = -Yk[i,j].imag
                    Ymnk[1,0] = +Yk[i,j].imag
                    Ymnk[1,1] = +Yk[i,j].real

                    # put it in the Y matrix
                    I = self.Kk * i + 2 * k
                    J = self.Kk * j + 2 * k
                    Y[I:I+2,J:J+2] = Ymnk

        return Y

    def calc_V0(self):

        # create node voltages using AC simulation results
        ac = AC('HB.AC', start=self.freqs[1], stop=self.freqs[self.K], numpts=self.K)
        ac.run(self.netlist)
        vac = ac.get_ac_solution()
        vdc = ac.get_dc_solution()

        V = np.zeros((self.W, 1))
        vt = np.zeros((self.N, self.S))
        for i in range(self.N):

            vf = np.zeros(self.K+1, dtype=complex)
            vf[0] = vdc[i]
            for k in range(1, self.K+1):
                vf[k] = vac[k-1,i]
                vf[k] = 0 # remove effect of AC for now

            for s in range(self.S):
                vt[i,s] = vf[0].real + 2 * (vf[1].real * np.cos(2. * np.pi * 1 * s / self.S) -
                                            vf[1].imag * np.sin(2. * np.pi * 1 * s / self.S))
                # vt[i,s] = min([vt[i,s], +0.4])
                # vt[i,s] = max([vt[i,s], -0.4])

            # plt.plot(vt[i,:])
            # plt.show()

            for k in range(self.K+1):
                vk = 0
                for s in range(self.S):
                    vk = vk + vt[i,s] * (np.cos(2. * np.pi * k * s / self.S) -
                                         1j * np.sin(2. * np.pi * k * s / self.S))
                vk = vk / self.S
                V[self.Kk*i+2*k] += vk.real
                V[self.Kk*i+2*k+1] += vk.imag

        return V, vt

    def ifft_V(self, V, vt):
        vt[:,:] = 0.
        for i in range(self.N):

            # assemble complex array of spectra for node 'i'
            vf = np.zeros(self.K+1, dtype=complex)
            for k in range(self.K+1):
                vf[k] = V[self.Kk*i+2*k+0] + 1j * V[self.Kk*i+2*k+1]

            # compute inverse fourier transform of voltage waveform
            for s in range(self.S):
                vt[i,s] = vf[0].real
                for k in range(1, self.K+1):
                    vt[i,s] = vt[i,s] + 2 * (vf[k].real * np.cos(2. * np.pi * k * s / self.S) -
                                             vf[k].imag * np.sin(2. * np.pi * k * s / self.S))

    def fft_G(self, gt, G):
        G[:,:,:] = 0
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.K+1):
                    for s in range(self.S):
                        G[i,j,k] = G[i,j,k] + gt[i,j,s] * (np.cos(2. * np.pi * k * s / self.S) -
                                                           1j * np.sin(2. * np.pi * k * s / self.S))
                    G[i,j,k] = G[i,j,k] / self.S

    def fft_Inl(self, it, Inl):
        Inl[:,:] = 0
        for i in range(self.N):
            for k in range(self.K+1):
                Inlk = complex()
                for s in range(self.S):
                    Inlk = Inlk + it[i,s] * (np.cos(2. * np.pi * k * s / self.S) -
                                           1j * np.sin(2. * np.pi * k * s / self.S))
                Inlk = Inlk / self.S
                Inl[self.Kk*i+2*k] = Inlk.real
                Inl[self.Kk*i+2*k+1] = Inlk.imag





