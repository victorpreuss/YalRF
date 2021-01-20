import sys
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

class MultiToneHarmonicBalance:

    def __init__(self, name, freq, numharmonics):
        self.name = name
        self.freq = freq
        self.numharmonics = numharmonics

        self.options = options.copy()

    def get_node_idx(self, node):
        return self.netlist.get_node_idx(node) - 1

    def plot_v(self, node):
        idx = self.get_node_idx(node)
        plt.title('Frequency-domain $V(j \omega)$')
        plt.stem(self.freqs, abs(self.Vf[idx,:]), use_line_collection=True, markerfmt='^')
        for f, v in zip(self.freqs, self.Vf[idx,:]):
            label = r'${:.3f} \angle {:.1f}^\circ $'.format(abs(v), np.degrees(np.angle(v)))
            plt.annotate(label, (f,abs(v)), textcoords="offset points", xytext=(0,10), ha='center')
        plt.xlabel('frequency [Hz]')
        plt.ylabel('magnitude [V]')
        plt.grid()
        plt.tight_layout()

    def run(self, netlist):

        self.netlist = netlist.copy()

        # get device data from netlist
        self.lin_devs = self.netlist.get_linear_devices()
        self.nonlin_devs = self.netlist.get_nonlinear_devices()

        # print option to make large matrices readable
        np.set_printoptions(precision=4, threshold=sys.maxsize, linewidth=160)

        """ Get setup variables """

        self.f1 = self.freq[0]
        self.f2 = self.freq[1]

        self.K1 = self.numharmonics[0]
        self.K2 = self.numharmonics[1]
        self.K = (2 * self.K2 + 1) * self.K1 + self.K2

        self.l0 = self.f1 / (2 * self.K2 + 1)

        self.freqs = np.zeros((self.K+1,1))
        self.mapfreqs = np.zeros((self.K+1,1))
        i = 0
        for k1 in range(self.K1+1):
            for k2 in range(-self.K2, self.K2+1):
                if k1 == 0 and k2 < 0:
                    continue
                k = (2 * self.K2 + 1) * k1 + k2
                self.freqs[i] = k1 * self.f1 + k2 * self.f2
                self.mapfreqs[i] = k * self.l0
                i += 1
                
        self.S = 2 * self.K + 1
        self.N = self.netlist.get_num_nodes() - 1
        self.W = self.S * self.N
        
        # print('Freqs = {}'.format(self.freqs))

        self.DFT = np.ones((self.S,self.S))
        for k in range(self.K):
            for s in range(self.S):
                self.DFT[2*k+1,s] = np.cos(2*np.pi*(k+1)*s/self.S);
                self.DFT[2*k+2,s] = - np.sin(2*np.pi*(k+1)*s/self.S);
        self.DFT = self.DFT / self.S;

        self.IDFT = np.ones((self.S,self.S))
        for s in range(self.S):
            for k in range(self.K):
                self.IDFT[s,2*k+1] = 2 * np.cos(2*np.pi*(k+1)*s/self.S)
                self.IDFT[s,2*k+2] = - 2 * np.sin(2*np.pi*(k+1)*s/self.S)

        """ Independent current sources (Is) """

        Is = self.calc_Is()

        # print('Is = {}'.format(Is))

        """ Transadmittance matrix Y(jw) """

        Y = self.calc_Y()

        # print('Y = {}'.format(Y))
        
        """ Initial voltage estimation for each node V(jw) and v(t) """

        V = self.calc_V0()
        vt = self.ifft(V)

        # print('V = {}'.format(V))
        # print('Vt = {}'.format(vt))

        """ Run Harmonic Balance Solver """

        # V, converged = self.hb_loop(Is, Y, V, vt)

        # if first attempt fails, try continuation method
        converged = False
        if converged == False:
            V = self.calc_V0()
            vt = self.ifft(V)
            Vprev = V.copy

            maxiter = 100
            alpha = 0.1
            converged = False
            itercnt = 0
            while (not converged) and (itercnt < maxiter):

                Iscpy = Is.copy()
                for i in range(self.N):
                    for k in range(1, self.K+1):
                        Iscpy[self.S*i+2*k-1] *= alpha
                        Iscpy[self.S*i+2*k] *= alpha

                print('Alpha = {}'.format(alpha))

                V, issolved = self.hb_loop(Iscpy, Y, V, vt)

                print('Converged: {}\n'.format(issolved))
                
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
                    vt = self.ifft(V)
                    alpha = alpha / 1.1
                    if alpha <= 0.001:
                        converged = False
                        break

                itercnt += 1

            print('End convergence = {}'.format(converged))
            print('Source Stepping iterations = {}'.format(itercnt))

            n = self.netlist.get_node_idx('n2') - 1
            print('{} {}'.format(np.abs(self.freqs[0]), V[n*self.S]))
            for k in range(1, self.K+1):
                i = n * self.S + 2 * (k - 1) + 1
                v = V[i] + 1j * V[i+1]
                print('{} {} {}'.format(np.abs(self.freqs[k]), np.abs(v), np.degrees(np.angle(v))))

        if not converged:
            return False, None, None, None, None

        return False, None, None, None, None
        # return converged, self.freqs, self.Vf, self.time, self.Vt

    def hb_loop(self, Is, Y, V, vt):

        it = np.zeros((self.W,1))
        dIdV = np.zeros((self.W,self.W))

        converged = False
        maxiter = self.options['maxiter']
        itercnt = 0
        while True:

            """ Time-domain v(t) """

            # vtold = vt.copy()

            # for k in range(1, self.K+1):
            #     if self.freqs[k] < 0:
            #         for n in range(self.N):
            #             i = n * self.S + 2 * (k - 1) + 2
            #             V[i] = - V[i]

            vt = self.ifft(V)

            # print('vt = {}'.format(vt))

            """ Time-domain g(t) and i(t) waveforms """

            # run a time-varying oppoint analysis on the nonlinear devices
            dIdV[:,:] = 0
            it[:] = 0
            for dev in self.nonlin_devs:
                if isinstance(dev, Diode):
                    n1 = dev.n1 - 1
                    n2 = dev.n2 - 1
                    Id = np.zeros((self.S,1))
                    gd = np.zeros((self.S))
                    for s in range(self.S):
                        v1 = vt[n1*self.S+s] if n1 >= 0 else 0
                        v2 = vt[n2*self.S+s] if n2 >= 0 else 0
                        vd = v1 - v2
                        Id[s] = dev.get_i(vd)
                        gd[s] = dev.get_g(vd)

                    gf = self.DFT @ np.diag(gd) @ self.IDFT

                    if n1 >= 0:
                        i = n1 * self.S
                        j = i + self.S
                        dIdV[i:j,i:j] += gf
                        it[i:j] += Id
                    if n2 >= 0:
                        i = n2 * self.S
                        j = i + self.S
                        dIdV[i:j,i:j] += gf
                        it[i:j] -= Id
                    if n1 >= 0 and n2 >= 0:
                        i = n1 * self.S
                        j = n2 * self.S
                        dIdV[i:i+self.S,j:j+self.S] -= gf
                        dIdV[j:j+self.S,i:i+self.S] -= gf

            """ Calculating the Jacobian Matrix J(jw) """

            J = Y + dIdV # + Omega * dQdV

            # print('J = {}'.format(J))

            """ Linear current Il(jw) """

            Il = Y @ V + Is

            """ Nonlinear current Inl(jw) """

            Inl = self.fft(it)

            # for k in range(1, self.K+1):
            #     if self.freqs[k] < 0:
            #         for n in range(self.N):
            #             i = n * self.S + 2 * (k - 1) + 2
            #             Inl[i] = - Inl[i]

            """ Calculate error function F(jw) """

            F = Il + Inl # + j * Omega * Q

            # print('F = {}'.format(F))

            """ Check algorithm termination conditions """

            converged = self.hb_converged(Il, Inl)
            itercnt += 1
            if (converged and itercnt >= 3) or itercnt >= maxiter:
                print('HB error = {:.2e}'.format(np.sum(np.abs(F))))
                print('Total number of iterations: {}'.format(itercnt))
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
            # V = V - scipy.linalg.inv(J) @ F

            # print('V = {}'.format(V))

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
        S = self.S
        Is = np.zeros((self.W, 1))
        for dev in self.lin_devs:
            if isinstance(dev, CurrentSource):
                iac = dev.ac * np.exp(1j * dev.phase)
                
                if dev.freq == self.f1:
                    fidx = 4 * self.K2 + 1
                else:
                    fidx = 1

                n1 = dev.n1 - 1
                if n1 >= 0:
                    Is[S*n1+0] += dev.dc
                    Is[S*n1+fidx] += iac.real / 2.
                    Is[S*n1+fidx+1] += iac.imag / 2.

                n2 = dev.n2 - 1
                if n2 >= 0:
                    Is[S*n2+0] -= dev.dc
                    Is[S*n2+fidx] -= iac.real / 2.
                    Is[S*n2+fidx+1] -= iac.imag / 2.

        return Is

    def calc_Y(self):
        Y = np.zeros((self.W, self.W))
        for dev in self.lin_devs:

            # add stamps for DC first
            if isinstance(dev, Resistor):
                n = (dev.n1 - 1) * self.S
                m = (dev.n2 - 1) * self.S
                y = 1. / dev.R
                if (dev.n1 > 0):
                    Y[n,n] += y
                if (dev.n2 > 0):
                    Y[m,m] += y
                if (dev.n1 > 0 and dev.n2 > 0):
                    Y[n,m] -= y
                    Y[m,n] -= y
            elif isinstance(dev, Gyrator):
                n1 = (dev.n1 - 1) * self.S
                n2 = (dev.n2 - 1) * self.S
                n3 = (dev.n3 - 1) * self.S
                n4 = (dev.n4 - 1) * self.S
                if (dev.n1 > 0 and dev.n2 > 0):
                    Y[n1,n2] += 1
                    Y[n2,n1] -= 1
                if (dev.n1 > 0 and dev.n3 > 0):
                    Y[n1,n3] -= 1
                    Y[n3,n1] += 1
                if (dev.n2 > 0 and dev.n4 > 0):
                    Y[n2,n4] += 1
                    Y[n4,n2] -= 1
                if (dev.n3 > 0 and dev.n4 > 0):
                    Y[n3,n4] -= 1
                    Y[n4,n3] += 1

            for k in range(1, self.K+1):
                if isinstance(dev, Resistor):
                    n = (dev.n1 - 1) * self.S + 2 * (k - 1) + 1
                    m = (dev.n2 - 1) * self.S + 2 * (k - 1) + 1
                    y = 1. / dev.R
                    Ymnk = np.array([[y, 0],[0, y]])
                    if (dev.n1 > 0):
                        Y[n:n+2,n:n+2] += Ymnk 
                    if (dev.n2 > 0):
                        Y[m:m+2,m:m+2] += Ymnk 
                    if (dev.n1 > 0 and dev.n2 > 0):
                        Y[n:n+2,m:m+2] -= Ymnk 
                        Y[m:m+2,n:n+2] -= Ymnk 
                elif isinstance(dev, Gyrator):
                    n1 = (dev.n1 - 1) * self.S + 2 * (k - 1) + 1
                    n2 = (dev.n2 - 1) * self.S + 2 * (k - 1) + 1
                    n3 = (dev.n3 - 1) * self.S + 2 * (k - 1) + 1
                    n4 = (dev.n4 - 1) * self.S + 2 * (k - 1) + 1
                    Ymnk = np.array([[1, 0],[0, 1]])
                    if (dev.n1 > 0 and dev.n2 > 0):
                        Y[n1:n1+2,n2:n2+2] += Ymnk
                        Y[n2:n2+2,n1:n1+2] -= Ymnk
                    if (dev.n1 > 0 and dev.n3 > 0):
                        Y[n1:n1+2,n3:n3+2] -= Ymnk
                        Y[n3:n3+2,n1:n1+2] += Ymnk
                    if (dev.n2 > 0 and dev.n4 > 0):
                        Y[n2:n2+2,n4:n4+2] += Ymnk
                        Y[n4:n4+2,n2:n2+2] -= Ymnk
                    if (dev.n3 > 0 and dev.n4 > 0):
                        Y[n3:n3+2,n4:n4+2] -= Ymnk
                        Y[n4:n4+2,n3:n3+2] += Ymnk

        return Y

    def calc_V0(self):
        # TODO: create better initial solution
        X = DC('HB.DC').run(self.netlist)
        V = np.zeros((self.W,1))

        return V

    def ifft(self, X):
        xt = np.zeros((self.W,1))
        for n in range(self.N):
            i = n * self.S
            j = i + self.S
            xt[i:j] = self.IDFT @ X[i:j]
        
        return xt

    def fft(self, xt):
        X = np.zeros((self.W,1))
        for n in range(self.N):
            i = n * self.S
            j = i + self.S
            X[i:j] = self.DFT @ xt[i:j]
        
        return X


