import numpy as np
from scipy import linalg

from yalrf.Netlist import Netlist
from yalrf.Analyses import AC, DC
from yalrf.Devices import *
from yalrf.Utils import hb_logger as logger

import logging
logger.setLevel(logging.INFO)

import matplotlib.pyplot as plt

class HarmonicBalance:

    def __init__(self, name, freq, numharmonics=4):
        self.name = name
        self.freq = freq
        self.numharmonics = numharmonics

    def run(self, netlist):

        # get data from main netlist
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

        S = 2 * K + 1                   # time array length for DFT/IDFT
        dt = T / S                      # timestep in the transient waveform

        # array with the Kth harmonic frequencies
        freqs = self.freq * np.linspace(0, K, K+1)

        """ Independent current sources (Is) """

        Is = np.zeros((W, 1))
        for dev in lin_devs:
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

        print('Is = {}'.format(Is))

        """ Transadmittance matrix (Y) """

        Y = np.zeros((W, W))
        for k in range(K+1):

            # add ground node
            Yk = np.zeros((N+1, N+1), dtype=complex)

            # add currently supported devices stamps
            for dev in lin_devs:
                if isinstance(dev, (Resistor, Capacitor, Inductor, Gyrator)):
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
        vt = np.zeros((N, S))
        for i in range(N):

            vf = np.zeros(K+1, dtype=complex)
            vf[0] = vdc[i]
            for k in range(1, K+1):
                vf[k] = vac[k-1,i]
                vf[k] = 0 # remove effect of AC for now

            for s in range(S):
                vt[i,s] = vf[0].real + 2 * (vf[1].real * np.cos(2. * np.pi * 1 * s / S) -
                                            vf[1].imag * np.sin(2. * np.pi * 1 * s / S))
                # vt[i,s] = min([vt[i,s], +0.4])
                # vt[i,s] = max([vt[i,s], -0.4])

            # plt.plot(vt[i,:])
            # plt.show()

            for k in range(K+1):
                vk = 0
                for s in range(S):
                    vk = vk + vt[i,s] * (np.cos(2. * np.pi * k * s / S) -
                                         1j * np.sin(2. * np.pi * k * s / S))
                vk = vk / S
                V[Kk*i+2*k] += vk.real
                V[Kk*i+2*k+1] += vk.imag

        print('V = {}'.format(V))

        for a in range(30): # run 30 iterations of HB

            """ Time-domain v(t) """

            vtold = vt.copy()
            vt = np.zeros((N, S))
            for i in range(N):

                # assemble complex array of spectra for node 'i'
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
            gt = np.zeros((N+1, N+1, S))
            it = np.zeros((N+1, S))
            for s in range(S):

                # get the conductances from DC oppoint
                for dev in nonlin_devs:
                    if isinstance(dev, Diode):

                        # get diode nodes
                        n1 = dev.n1
                        n2 = dev.n2

                        # need to get the correct previous voltage for limiting
                        # TODO: unuderstande if this is needed
                        v1old = vtold[n1-1,s] if n1 > 0 else 0
                        v2old = vtold[n2-1,s] if n2 > 0 else 0
                        Vdold = v1old - v2old

                        # calculate current over the diode
                        v1 = vt[n1-1,s] if n1 > 0 else 0
                        v2 = vt[n2-1,s] if n2 > 0 else 0
                        Vd = v1 - v2

                        Id = dev.get_i(Vd, Vdold)
                        gd = dev.get_g(Vd, Vdold)

                        # USING SECANT METHOD
                        # Id2 = dev.get_i(Vd * 1.01, Vdold)
                        # gd = (Id2 - Id) / (0.01 * Vd)

                        gt[n1,n1,s] = gt[n1,n1,s] + gd
                        gt[n2,n2,s] = gt[n2,n2,s] + gd
                        gt[n1,n2,s] = gt[n1,n2,s] - gd
                        gt[n2,n1,s] = gt[n2,n1,s] - gd

                        it[n1,s] = it[n1,s] + Id
                        it[n2,s] = it[n2,s] - Id

                        # USING SECANT METHOD
                        # # apply perturbation to node n1
                        # if n1 >= 0:
                        #     Vd1 = 1.01 * v1 - v2
                        #     Id1 = dev.get_i(Vd1, Vdold)
                        #     gd1 = (Id1 - Id) / (0.01 * v1)

                        #     gt[n1,n1,s] = gt[n1,n1,s] + gd1
                        #     it[n1,s] = it[n1,s] + Id

                        #     if n2 >= 0:
                        #         gt[n2,n1,s] = gt[n2,n1,s] - gd1

                        # # apply perturbation to node n2
                        # if n2 >= 0:
                        #     Vd2 = v1 - 1.01 * v2
                        #     Id2 = dev.get_i(Vd2, Vdold)
                        #     gd2 = (Id2 - Id) / (- 0.01 * v2)

                        #     gt[n2,n2,s] = gt[n2,n2,s] + gd2
                        #     it[n2,s] = it[n2,s] - Id

                        #     if n1 >= 0:
                        #         gt[n1,n2,s] = gt[n1,n2,s] - gd2

                    if isinstance(dev, BJT):

                        # get BJT nodes
                        B = dev.n1 
                        C = dev.n2 
                        E = dev.n3 

                        # calculate currents at the BJT
                        Vb = vt[B-1,s] if B > 0 else 0
                        Vc = vt[C-1,s] if C > 0 else 0
                        Ve = vt[E-1,s] if E > 0 else 0
                        Vs = 0

                        Ib, Ic, Ie = dev.get_i(Vb, Vc, Ve, Vs)
                        gmu, gpi, gmf, gmr = dev.get_g(Vb, Vc, Ve, Vs)

                        gt[B,B,s] = gt[B,B,s] + gmu + gpi
                        gt[B,C,s] = gt[B,C,s] - gmu
                        gt[B,E,s] = gt[B,E,s] - gpi
                        gt[C,B,s] = gt[C,B,s] - gmu + gmf + gmr
                        gt[C,C,s] = gt[C,C,s] + gmu - gmr
                        gt[C,E,s] = gt[C,E,s] - gmf
                        gt[E,B,s] = gt[E,B,s] - gpi - gmf - gmr
                        gt[E,C,s] = gt[E,C,s] + gmr
                        gt[E,E,s] = gt[E,E,s] + gpi + gmf

                        it[B,s] = it[B,s] + Ib
                        it[C,s] = it[C,s] + Ic
                        it[E,s] = it[E,s] - Ie

                        # USING SECANT METHOD (need to understand better what to do here)
                        # Ibn, Icn, Ien = dev.get_i(Vb * 1.01, Vc, Ve, Vs)

                        # gbb = (Ibn - Ib) / (0.01 * Vb)
                        # gcb = (Icn - Ic) / (0.01 * Vb)
                        # geb = (Ien - Ie) / (0.01 * Vb)

                        # gt[B,B,s] = gt[B,B,s] + gbb
                        # gt[C,B,s] = gt[C,B,s] + gcb
                        # gt[E,B,s] = gt[E,B,s] + geb

                        # Ibn, Icn, Ien = dev.get_i(Vb, Vc * 1.01, Ve, Vs)

                        # gcc = (Icn - Ic) / (0.01 * Vc)
                        # gbc = (Ibn - Ib) / (0.01 * Vc)
                        # gec = (Ien - Ie) / (0.01 * Vc)

                        # gt[C,C,s] = gt[C,C,s] + gcc
                        # gt[B,C,s] = gt[B,C,s] + gbc
                        # gt[E,C,s] = gt[E,C,s] + gec

                        # Ibn, Icn, Ien = dev.get_i(Vb, Vc, Ve * 1.01, Vs)

                        # gee = (Ien - Ie) / (0.01 * Ve)
                        # gbe = (Ibn - Ib) / (0.01 * Ve)
                        # gce = (Icn - Ic) / (0.01 * Ve)

                        # gt[E,E,s] = gt[E,E,s] + gee
                        # gt[B,E,s] = gt[B,E,s] + gbe
                        # gt[C,E,s] = gt[C,E,s] + gce

                        # it[B,s] = it[B,s] + Ib
                        # it[C,s] = it[C,s] + Ic
                        # it[E,s] = it[E,s] - Ie

            # remove ground node
            gt = gt[1:,1:,:]
            it = it[1:,:]

            # print('vt = {}'.format(vt))
            # print('gt = {}'.format(gt))
            # print('it = {}'.format(it))

            # compute the spectrum of each port conductance
            G = np.zeros((N, N, K+1), dtype=complex)
            for i in range(N):
                for j in range(N):
                    for k in range(K+1):
                        for s in range(S):
                            G[i,j,k] = G[i,j,k] + gt[i,j,s] * (np.cos(2. * np.pi * k * s / S) -
                                                               1j * np.sin(2. * np.pi * k * s / S))
                        G[i,j,k] = G[i,j,k] / S

            """ Nonlinear current (Ig) """

            # TODO: generalize this section for 3-port devices
            Ig = np.zeros((W, 1))
            for i in range(N):
                for k in range(K+1):
                    Igk = complex()
                    for s in range(S):
                        Igk = Igk + it[i,s] * (np.cos(2. * np.pi * k * s / S) -
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
            J = Y + dIdV # + Omega * dQdV

            for i in range(N):
                k = Kk * i
                J[k+1,k+1] = 1.

            F = Y @ V + Ig - Is # + j * Omega * Q
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



