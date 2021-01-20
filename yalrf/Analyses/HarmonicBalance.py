import numpy as np
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

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

    def __init__(self, name, freq, numharmonics=5):
        self.name = name
        self.freq = freq
        self.numharmonics = numharmonics

        self.options = options.copy()

    def get_node_idx(self, node):
        return self.netlist.get_node_idx(node) - 1

    def plot_v(self, node):
        idx = self.get_node_idx(node)
        plt.figure(figsize=(10,8))

        plt.subplot(211)
        plt.title('Time-domain $V(t)$')
        plt.plot(self.time, self.Vt[idx,:])
        plt.xlabel('time [s]')
        plt.ylabel('V' + node + ' [V]')
        plt.grid()

        plt.subplot(212)
        plt.title('Frequency-domain $V(j \omega)$')
        plt.stem(self.freqs, abs(self.Vf[idx,:]), use_line_collection=True, markerfmt='^')
        for f, v in zip(self.freqs, self.Vf[idx,:]):
            label = r'${:.3f} \angle {:.1f}^\circ $'.format(abs(v), np.degrees(np.angle(v)))
            plt.annotate(label, (f,abs(v)), textcoords="offset points", xytext=(0,10), ha='center')
        plt.xlabel('frequency [Hz]')
        plt.ylabel('magnitude [V]')
        plt.grid()
        plt.tight_layout()

    def run(self, netlist, V0=None):

        self.netlist = netlist.copy()

        # get device data from netlist
        self.lin_devs = self.netlist.get_linear_devices()
        self.nonlin_devs = self.netlist.get_nonlinear_devices()

        # print option to make large matrices readable
        np.set_printoptions(precision=4, suppress=False, linewidth=150)

        """ Get basic sizes """

        self.K = self.numharmonics           # frequency of the 1st harmonic
        self.N = netlist.get_num_nodes() - 1 # total number of nodes in the netlist without 'gnd'

        self.Kk = 2 * (self.K + 1)           # size of block matrices in the system
        self.W = self.Kk * self.N            # total size of the matrix system 2 * (K+1) * N
        
        self.T = 1. / self.freq              # period of the 1st harmonic
        self.S = 2 * self.K + 1              # time array length for DFT/IDFT

        # array with the Kth harmonic frequencies
        self.freqs = self.freq * np.linspace(0, self.K, self.K+1)

        """ Independent current sources (Is) """

        Is = self.calc_Is()

        # print('Is = {}'.format(Is))

        """ Transadmittance matrix Y(jw) """

        Y = self.calc_Y()

        # print('Y = {}'.format(Y))
        
        """ Initial voltage estimation for each node V(jw) and v(t) """

        if V0 is None:
            V, vt = self.calc_V0()
        else:
            V = V0
            vt = np.zeros((self.N, self.S))
            self.ifft_V(V, vt)

        # print('V = {}'.format(V))
        # print('Vt = {}'.format(vt))

        """ Run Harmonic Balance Solver """

        # V, converged = self.hb_loop(Is, Y, V, vt)
        # print(converged)
        converged = False

        # if first attempt fails, try continuation method
        if converged == False:
            # V, vt = self.calc_V0()
            Vprev = V.copy()

            if V0 is None:
                alpha = 0.1
            else:
                alpha = 1.
            maxiter = 100
            converged = False
            itercnt = 0
            while (not converged) and (itercnt < maxiter):

                Iscpy = Is.copy()
                for i in range(self.N):
                    for k in range(1, self.K+1):
                        Iscpy[self.Kk*i+2*k+0] *= alpha
                        Iscpy[self.Kk*i+2*k+1] *= alpha

                V, issolved = self.hb_loop(Iscpy, Y, V, vt)

                print('alpha = {}'.format(alpha))
                print('issolved = {}'.format(issolved))
                
                if issolved:
                    
                    if alpha >= 1.0:
                        converged = True
                        break

                    Vprev = V.copy()
                    alpha = 1.25 * alpha
                    if alpha >= 1.0:
                        alpha = 1.0

                else:
                    V = Vprev
                    self.ifft_V(V, vt)
                    alpha = alpha / 1.1
                    if alpha <= 0.001:
                        converged = False
                        break

                itercnt += 1

            print('converged = {}'.format(converged))
            print('iterations = {}'.format(itercnt))

        if not converged:
            return converged, None, None, None, None

        S = 32 * self.K # increase number of time samples
        Vt = np.zeros((self.N, S))
        Vf = np.zeros((self.N, self.K+1), dtype=complex)
        for i in range(self.N):

            # assemble complex array of spectra for node 'i'
            for k in range(self.K+1):
                Vf[i,k] = V[self.Kk*i+2*k+0] + 1j * V[self.Kk*i+2*k+1]

            # compute inverse fourier transform of voltage waveform
            for s in range(S):
                Vt[i,s] = Vf[i,0].real
                for k in range(1, self.K+1):
                    Vt[i,s] = Vt[i,s] + 2 * (Vf[i,k].real * np.cos(2. * np.pi * k * s / (S / 2)) -
                                             Vf[i,k].imag * np.sin(2. * np.pi * k * s / (S / 2)))


        self.time = 2 * self.T / S * np.linspace(0, S, S)
        self.Vt = Vt
        self.Vf = Vf

        # store answer array
        self.X = V

        return converged, self.freqs, self.Vf, self.time, self.Vt

    def hb_loop(self, Is, Y, V, vt):
        gt = np.zeros((self.N+1, self.N+1, self.S))             # time-varying nonlinear conductances
        G = np.zeros((self.N, self.N, self.K+1), dtype=complex) # frequency-domain nonlinear conductances

        it = np.zeros((self.N+1, self.S)) # time-domain nonlinear current
        Inl = np.zeros((self.W, 1))       # frequency-domain nonlinear current

        # matrices that form dI/dV
        Toe = np.zeros((self.W, self.W))
        Han = np.zeros((self.W, self.W))
        D = np.eye(self.W)

        # DC elements are 1/2 (from DFT derivation)
        for i in range(self.N):
            k = self.Kk * i
            D[k,k] = 0.5
            D[k+1,k+1] = 0.5

        converged = False
        maxiter = self.options['maxiter']
        itercnt = 0
        while True:

            """ Time-domain v(t) """

            # vtold = vt.copy()
            self.ifft_V(V, vt)

            # print('vt = {}'.format(vt))

            """ Time-domain g(t) and i(t) waveforms """

            # run a time-varying oppoint analysis on the nonlinear devices
            gt[:,:,:] = 0
            it[:,:] = 0
            for s in range(self.S):
                for dev in self.nonlin_devs:
                    dev.add_hb_stamps(vt, it, gt, s)

            # print('gt = {}'.format(gt))
            # print('it = {}'.format(it))
            
            """ Nonlinear conductance G(jw) """

            self.fft_G(gt[1:,1:,:], G)

            # print('G = {}'.format(G))

            """ Creating the dIdV matrix from G(jw) """

            dIdV = self.calc_dIdV(G, Toe, Han, D)

            # print('dIdV = {}'.format(dIdV))

            """ Calculating the Jacobian Matrix J(jw) """

            J = Y + dIdV # + Omega * dQdV

            # need to add a dummy value to the complex part of the DC voltages
            for i in range(self.N):
                k = self.Kk * i
                J[k+1,k+1] = 1.

            # print('J = {}'.format(J))

            """ Linear current Il(jw) """

            Il = Y @ V - Is

            """ Nonlinear current Inl(jw) """

            self.fft_Inl(it[1:,:], Inl)

            """ Calculate error function F(jw) """

            F = Il + Inl # + j * Omega * Q

            # print('F = {}'.format(F))

            """ Check algorithm termination conditions """

            converged = self.hb_converged(Il, Inl)
            itercnt += 1
            if (converged and itercnt >= 3) or itercnt >= maxiter:
                # print('F = {}'.format(F))
                return V, converged

            """ Calculate next voltage guess using NR """

            # currently just using standard LU factorization
            if False:
                Jn = scipy.sparse.csc_matrix(J)
                lu = scipy.sparse.linalg.splu(Jn)
                dV = lu.solve(F)
            else:
                lu, piv = scipy.linalg.lu_factor(J)
                dV = scipy.linalg.lu_solve((lu, piv), F)

            V = V - dV
            # V = V - linalg.inv(J) @ F

            # ensure there is no complex part on the DC voltages
            for i in range(self.N):
                V[self.Kk*i+1] = 0.

            # print('V = {}'.format(V))

        print('Total number of iterations: {}'.format(itercnt))

    def hb_converged(self, Il, Inl):
        abstol = self.options['abstol']
        reltol = self.options['reltol']
        converged = True

        F = Il + Inl
        Frel = 2 * abs(F / (Il - Inl + 1e-15))
        for i in range(self.W):
            if abs(F[i]) > abstol and Frel[i] > reltol:
                converged = False
                break

        return converged

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
                    Is[Kk*n1+3] -= iac.imag / 2.

                n2 = dev.n2 - 1
                if n2 >= 0:
                    Is[Kk*n2+0] -= dev.dc
                    Is[Kk*n2+1] -= 0.
                    Is[Kk*n2+2] -= iac.real / 2.
                    Is[Kk*n2+3] += iac.imag / 2.

        return Is

    def calc_Y(self):
        Y = np.zeros((self.W, self.W))
        for k in range(self.K+1):

            # add ground node
            Yk = np.zeros((self.N+1, self.N+1), dtype=complex)

            # add currently supported devices stamps
            for dev in self.lin_devs:
                if isinstance(dev, (Resistor, Capacitor, Inductor, Gyrator, IdealHarmonicFilter)):
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

        # plt.plot(vt[0,:])
        # plt.show()

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

    def calc_dIdV(self, G, Toe, Han, D):
        N = self.N
        K = self.K
        Kk = self.Kk

        # populating Toeplitz and Hankel matrices
        Toe[:,:] = 0
        Han[:,:] = 0
        for i in range(N):
            for j in range(N):
                Tmn = np.zeros((Kk, Kk))
                Hmn = np.zeros((Kk, Kk))
                for k in range(K+1):
                    for l in range(k, K+1):
                        
                        # Toeplitz submatrix
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

                    for l in range(K+1):

                        # Hankel submatrix
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

                # insert Tmn and Hmn submatrices into the proper position
                I = Kk * i
                J = Kk * j
                Toe[I:I+Kk,J:J+Kk] = Tmn
                Han[I:I+Kk,J:J+Kk] = Hmn

        return (Toe + Han) @ D

