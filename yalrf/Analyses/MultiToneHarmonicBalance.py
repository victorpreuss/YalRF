import sys
import numpy as np
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
from scipy import optimize

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

    def __init__(self, name, freq=1e9, numharmonics=10):
        self.name = name
        self.freq = freq
        self.numharmonics = numharmonics

        self.options = options.copy()

    def get_node_idx(self, node):
        return self.netlist.get_node_idx(node) - 1

    def plot_v(self, node):
        n = self.get_node_idx(node)

        if self.freqs[1] > 1e9:
            funit = 'GHz'
            c = 1e9
        elif self.freqs[1] > 1e6:
            funit = 'MHz'
            c = 1e6
        elif self.freqs[1] > 1e3:
            funit = 'kHz'
            c = 1e3
        else:
            funit = 'Hz'
            c = 1

        plt.rc('axes', titlesize=14)    # fontsize of the axes title
        plt.rc('axes', labelsize=12)    # fontsize of the x and y labels

        if self.num_tones == 1:
            plt.figure(figsize=(14,5))
            plt.subplot(121)
        else:
            plt.figure(figsize=(10,6))

        # plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        # plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        # plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        # plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        # plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

        plt.title('Frequency-domain $V(j \omega)$')
        plt.stem(self.freqs / c, np.abs(self.Vf[n]), use_line_collection=True, markerfmt='r^', linefmt='r')
        for f, v in zip(self.freqs, self.Vf[n]):
            label = r'{:.3f} $\angle$ {:.1f}$^\circ$'.format(np.abs(v), np.degrees(np.angle(v)))
            # if f <= 3 * self.freqs[1]:
            plt.annotate(label, (f / c, np.abs(v)), textcoords="offset points", xytext=(0,5), ha='left', va='bottom', rotation=45)
        plt.xlabel('frequency [' + funit + ']')
        plt.ylabel('V(\'' + node + '\') [V]')
        plt.grid()
        plt.tight_layout()

        # if single-tone analysis, also plot time-domain waveform
        if self.num_tones == 1:
            S = 32 * self.K # increase number of time samples
            time = 2 / self.freq / S * np.linspace(0, S, S)
            Vt = np.zeros(S)
            for s in range(S):
                Vt[s] = self.Vf[n,0].real
                for k in range(1, self.K+1):
                    Vt[s] = Vt[s] + self.Vf[n,k].real * np.cos(2. * np.pi * k * s / (S / 2)) - \
                                    self.Vf[n,k].imag * np.sin(2. * np.pi * k * s / (S / 2))

            if time[-1] > 1:
                tunit = 's'
                c = 1
            elif time[-1] > 1e-3:
                tunit = 'ms'
                c = 1e3
            elif time[-1] > 1e-6:
                tunit = '$\mu$s'
                c = 1e6
            elif time[-1] > 1e-9:
                tunit = 'ns'
                c = 1e9
            else: # time[-1] > 1e-12
                tunit = 'ps'
                c = 1e12

            plt.subplot(122)
            plt.title('Time-domain $V(t)$')
            plt.plot(time * c, Vt)
            plt.xlabel('time [' + tunit + ']')
            plt.ylabel('V(\'' + node + '\') [V]')
            plt.grid()
            plt.tight_layout()

    def print_v(self, node):
        n = self.get_node_idx(node)
        print('Voltage at node: ' + node)
        print('Freq [Hz]\tVmag [V]\tPhase [°]')
        for k in range(0, self.K+1):
            v = self.Vf[n,k]
            print('{:.2e}\t{:.3e}\t{:.1f}'.format(self.freqs[k], np.abs(v), np.degrees(np.angle(v))))

    def get_v(self, node):
        n = self.get_node_idx(node)
        v = self.Vf[n]
        return v

    def run_oscillator(self, netlist, f0, numharmonics, V0, node):

        netlist = netlist.copy()

        # config harmonic balance
        self.freq = f0
        self.numharmonics = numharmonics

        # vin oscprobe
        Voscprobe = netlist.add_iac(self.name + '.I.Oscprobe', self.name + '.nosc1', 'gnd', ac=V0, freq=f0)
        netlist.add_gyrator(self.name + '.G.Oscprobe', self.name + '.nosc1', self.name + '.nosc2', 'gnd', 'gnd', 1)

        # ideal harmonic filter of oscprobe
        Zoscprobe = netlist.add_idealharmonicfilter(self.name + 'IHF.Oscprobe', self.name + '.nosc2', node, f0)

        def objFunc(x, info):

            # get new solution candidates
            fosc = x[0]
            Vosc = x[1]

            # update oscprobe values
            Voscprobe.ac = Vosc
            Voscprobe.freq = fosc
            Zoscprobe.freq = fosc

            # run harmonic balance
            self.freq = fosc
            if info['itercnt'] > 0:
                converged, freqs, Vf, time, Vt = self.run(netlist, self.V)
            else:
                converged, freqs, Vf, time, Vt = self.run(netlist)

            # if HB failed to converge, return a bad convergence value to minimizer
            if not converged:
                return 1e2

            # get nodes of IdealHarmonicFilter to calculate current across it
            node = info['osc_node']
            n1 = self.get_node_idx(node)
            n2 = self.get_node_idx(self.name + '.nosc2')

            # mag(Yosc) is the objective function to be minimized
            Voscx = Vf[n1,1]                             # voltage across oscprobe
            Iosc  = (Vf[n1,1] - Vf[n2,1]) * Zoscprobe.g  # current through oscprobe
            Yosc  = Iosc / Voscx

            info['itercnt'] += 1
            print('\nIter\tFreq [Hz]\tVosc [V]\tmag(Yosc)')
            print('{}\t{:.8f}\t{:.8f}\t{:.2e}\n'.format(info['itercnt'], fosc, Vosc, abs(Yosc)))

            if abs(Yosc) < 1e-12:
                print("Optimization terminated succesfully! |Y| < 1e-12")
                info['result'] = x
                raise SmallEnoughException()

            return abs(Yosc)

        class SmallEnoughException(Exception):
            pass

        try:
            x0 = [f0, V0]
            args = ({'itercnt' : 0, 'osc_node' : node, 'result' : 0},)
            xopt = optimize.fmin(func     = objFunc,
                                 x0       = x0,
                                 args     = args,
                                 xtol     = 1e-3,
                                 ftol     = 1e-3,
                                 maxfun   = 200,
                                 disp     = True,
                                 retall   = False,
                                 full_output = True)[0]
        except SmallEnoughException:
            xopt = args[0]['result']

        # get solution 
        fosc = xopt[0]
        Vosc = xopt[1]

        # update oscprobe values
        Voscprobe.ac = Vosc
        Voscprobe.freq = fosc
        Zoscprobe.freq = fosc

        # run harmonic balance
        self.freq = fosc
        converged, freqs, Vf, time, Vt = self.run(netlist, self.V)

        print('Frequency of oscillation = {} Hz'.format(fosc))
        print('Oscillation amplitude = {} V'.format(Vosc))

        return converged, freqs, Vf, time, Vt

    def run(self, netlist, V0=None):

        self.netlist = netlist.copy()

        # get device data from netlist
        self.lin_devs = self.netlist.get_linear_devices()
        self.nonlin_devs = self.netlist.get_nonlinear_devices()

        # print option to make large matrices readable
        np.set_printoptions(precision=4, threshold=sys.maxsize, linewidth=160)

        """ Get setup variables and create AFM frequency array """

        if isinstance(self.freq, float):
            self.num_tones = 1
            self.f1 = self.freq
            self.K = self.numharmonics
            self.freqs = self.f1 * np.linspace(0, self.K, self.K+1)
        elif isinstance(self.freq, list) and len(self.freq) == 2:
            self.num_tones = 2
            self.f1 = self.freq[0]
            self.f2 = self.freq[1]

            self.K1 = self.numharmonics[0]
            self.K2 = self.numharmonics[1]
            self.K = (2 * self.K2 + 1) * self.K1 + self.K2

            self.l0 = self.f1 / (2 * self.K2 + 1)

            self.freqs = np.zeros((self.K+1))
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
        else:
            print('ERROR: only one and two-tones are currently supported')
                
        self.S = 2 * self.K + 1
        self.N = self.netlist.get_num_nodes() - 1
        self.W = self.S * self.N
        
        # print('Freqs = {}'.format(self.freqs))

        """ Create DFT and IDFT transformation matrices """

        # self.DFT = np.ones((self.S,self.S))
        # for k in range(self.K):
        #     for s in range(self.S):
        #         # self.DFT[2*k+1,s] = np.cos(2*np.pi*(k+1)*s/self.S);
        #         # self.DFT[2*k+2,s] = - np.sin(2*np.pi*(k+1)*s/self.S);
        #         self.DFT[2*k+1,s] = 2 * np.cos(2*np.pi*(k+1)*s/self.S);
        #         self.DFT[2*k+2,s] = 2 * np.sin(2*np.pi*(k+1)*s/self.S);
        # self.DFT = self.DFT / self.S;

        self.IDFT = np.ones((self.S,self.S))
        for s in range(self.S):
            for k in range(self.K):
                # self.IDFT[s,2*k+1] = 2 * np.cos(2*np.pi*(k+1)*s/self.S)
                # self.IDFT[s,2*k+2] = - 2 * np.sin(2*np.pi*(k+1)*s/self.S)
                self.IDFT[s,2*k+1] = np.cos(2*np.pi*(k+1)*s/self.S)
                self.IDFT[s,2*k+2] = np.sin(2*np.pi*(k+1)*s/self.S)

        self.DFT = np.linalg.inv(self.IDFT)

        # check if DFT matrices are consistent
        # print(np.max(self.DFT @ self.IDFT - np.eye(S)))

        """ Independent current sources (Is) """

        Is = self.calc_Is()

        # print('Is = {}'.format(Is))

        """ Transadmittance matrix Y(jw) """

        Y = self.calc_Y()

        # print('Y = {}'.format(Y))
        
        """ Initial voltage estimation for each node V(jw) and v(t) """

        if V0 is None:
            V = self.calc_V0()
        else:
            V = V0

        # print('V = {}'.format(V))

        """ Run Harmonic Balance Solver """

        if V0 is None:
            converged = False
        else:
            V, converged = self.hb_loop(Is, Y, V)
            if not converged:
                V = self.calc_V0()

        # if first attempt fails, try continuation method
        if converged == False:
            Vprev = V.copy()
            maxiter = 40
            alpha = 0.1
            inc = 1.4
            dec = 1.2
            converged = False
            itercnt = 0
            while (not converged) and (itercnt < maxiter):

                Isalpha = Is.copy()
                for n in range(self.N):
                    for k in range(1, self.K+1):
                        Isalpha[self.S*n+2*k-1] *= alpha
                        Isalpha[self.S*n+2*k] *= alpha

                print('Alpha level: {}'.format(alpha))

                V, iterconverged = self.hb_loop(Isalpha, Y, V)

                print('HB iteration converged: {}\n'.format(iterconverged))
                
                if iterconverged:
                    if alpha >= 1.0:
                        converged = True
                        break

                    Vprev = V.copy()
                    alpha = inc * alpha
                    if alpha >= 1.0:
                        alpha = 1.0

                else:
                    V = Vprev.copy()
                    vt = self.ifft(V)
                    if alpha >= 1.0:
                        dec = 1 + 1.5 * (dec - 1)
                        inc = 1 + (inc - 1) / 1.5
                    alpha = alpha / dec
                    if alpha <= 0.01:
                        converged = False
                        break

                itercnt += 1

            print('Final convergence: {}'.format(converged))
            print('Number of source stepping iterations: {}\n'.format(itercnt))

        if not converged:
            return False, None, None, None, None

        if self.num_tones == 2:
            self.freqs = np.abs(self.freqs)

        # keep backup of solution (it can be used as I.C. for another run)
        self.V = V

        # Vf is an array with the voltage phasor at each node for all frequencies
        self.Vf = np.zeros((self.N,self.K+1), dtype=complex)
        for n in range(self.N):
            self.Vf[n,0] = self.V[n*self.S,0]
            for k in range(1, self.K+1):
                i = n * self.S + 2 * (k - 1) + 1
                An  = np.sqrt(self.V[i,0]**2 + self.V[i+1,0]**2)
                phi = np.arctan2(self.V[i+1,0], self.V[i,0])
                self.Vf[n,k] = An * np.exp(1j * phi)

        return converged, self.freqs, self.Vf, None, None

    def hb_loop(self, Is, Y, V):

        vt = self.ifft(V)
        it = np.zeros((self.W,1))
        qt = np.zeros((self.W,1))
        dIdV = np.zeros((self.W,self.W))
        dQdV = np.zeros((self.W,self.W))
        Omega = np.zeros((self.W,self.W))

        for n in range(self.N):
            for k in range(1,self.K+1):
                i = n * self.S + 2 * (k - 1) + 1
                omega = 2 * np.pi * np.abs(self.freqs[k])
                y = np.array([[0, -omega],[omega, 0]])
                Omega[i:i+2,i:i+2] = y

        converged = False
        maxiter = self.options['maxiter']
        itercnt = 0
        while True:

            """ Time-domain v(t) """

            vtprev = vt.copy()
            vt = self.ifft(V)

            # print('vt = {}'.format(vt))

            """ Time-domain g(t) and i(t) waveforms """

            # run a time-varying oppoint analysis on the nonlinear devices
            dIdV[:,:] = 0
            it[:] = 0
            qt[:] = 0
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

                        v1 = vtprev[n1*self.S+s] if n1 >= 0 else 0
                        v2 = vtprev[n2*self.S+s] if n2 >= 0 else 0
                        vdold = v1 - v2
                        
                        gd[s], Id[s] = dev.get_mthb_params(vd, vdold)

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

                elif isinstance(dev, BJT):
                    B = dev.n1 - 1
                    C = dev.n2 - 1
                    E = dev.n3 - 1

                    Ib = np.zeros((self.S,1))
                    Ic = np.zeros((self.S,1))
                    Ie = np.zeros((self.S,1))

                    Qbe = np.zeros((self.S,1))
                    Qbc = np.zeros((self.S,1))
                    Qsc = np.zeros((self.S,1))

                    gmu = np.zeros((self.S))
                    gpi = np.zeros((self.S))
                    gmf = np.zeros((self.S))
                    gmr = np.zeros((self.S))

                    Cbe = np.zeros((self.S))
                    Cbc = np.zeros((self.S))
                    Cbebc = np.zeros((self.S))
                    Csc = np.zeros((self.S))

                    for s in range(self.S):
                        Vb = vt[B*self.S+s] if B >= 0 else 0
                        Vc = vt[C*self.S+s] if C >= 0 else 0
                        Ve = vt[E*self.S+s] if E >= 0 else 0
                        Vs = 0

                        Vbold = vtprev[B*self.S+s] if B >= 0 else 0
                        Vcold = vtprev[C*self.S+s] if C >= 0 else 0
                        Veold = vtprev[E*self.S+s] if E >= 0 else 0

                        Ib[s], Ic[s], Ie[s], Qbe[s], Qbc[s], Qsc[s], gmu[s], gpi[s], gmf[s], gmr[s], Cbc[s], Cbe[s], Cbebc[s], Csc[s] = dev.get_hb_params(Vb, Vc, Ve, Vs, s, Vbold, Vcold, Veold)

                    gmuf = self.DFT @ np.diag(gmu) @ self.IDFT
                    gpif = self.DFT @ np.diag(gpi) @ self.IDFT
                    gmff = self.DFT @ np.diag(gmf) @ self.IDFT
                    gmrf = self.DFT @ np.diag(gmr) @ self.IDFT

                    Cbef = self.DFT @ np.diag(Cbe) @ self.IDFT
                    Cbcf = self.DFT @ np.diag(Cbc) @ self.IDFT
                    Cbebcf = self.DFT @ np.diag(Cbebc) @ self.IDFT
                    Cscf = self.DFT @ np.diag(Csc) @ self.IDFT

                    if B >= 0:
                        i = B * self.S
                        j = i + self.S
                        dIdV[i:j,i:j] += (gmuf + gpif)
                        dQdV[i:j,i:j] += (Cbcf + Cbef + Cbebcf)
                        it[i:j] += Ib
                        qt[i:j] += (Qbe + Qbc)

                        if C >= 0:
                            x = C * self.S
                            y = x + self.S
                            dIdV[i:j,x:y] += (- gmuf)
                            dIdV[x:y,i:j] += (- gmuf + gmff + gmrf)
                            dQdV[i:j,x:y] += (- Cbcf - Cbebcf)
                            dQdV[x:y,i:j] += (- Cbcf)

                        if E >= 0:
                            x = E * self.S
                            y = x + self.S
                            dIdV[i:j,x:y] += (- gpif)
                            dIdV[x:y,i:j] += (- gpif - gmff - gmrf)
                            dQdV[i:j,x:y] += (- Cbef)
                            dQdV[x:y,i:j] += (- Cbef - Cbebcf)

                    if C >= 0:
                        i = C * self.S
                        j = i + self.S
                        dIdV[i:j,i:j] += (+ gmuf - gmrf)
                        dQdV[i:j,i:j] += (+ Cbcf + Cscf)
                        it[i:j] += Ic
                        qt[i:j] -= Qbc - Qsc

                        if E >= 0:
                            x = E * self.S
                            y = x + self.S
                            dIdV[i:j,x:y] += (- gmff)
                            dIdV[x:y,i:j] += (+ gmrf)
                            dQdV[i:j,x:y] += 0.
                            dQdV[x:y,i:j] += (+ Cbebcf)

                    if E >= 0:
                        i = E * self.S
                        j = i + self.S
                        dIdV[i:j,i:j] += (+ gpif + gmff)
                        dQdV[i:j,i:j] += (+ Cbef)
                        it[i:j] -= Ie
                        qt[i:j] -= Qbe

            """ Calculating the Jacobian Matrix J(jw) """

            J = Y + dIdV + Omega @ dQdV

            # print('J = {}'.format(J))

            """ Linear current Il(jw) """

            Il = Y @ V - Is

            """ Nonlinear current Inl(jw) """

            Inl = self.fft(it) + Omega @ self.fft(qt)

            """ Calculate error function F(jw) """

            F = Il + Inl

            # print('F = {}'.format(F))

            """ Check algorithm termination conditions """

            converged = self.hb_converged(Il, Inl)
            itercnt += 1
            if (converged and itercnt >= 3) or itercnt >= maxiter:
                print('HB total error: {:.2e}'.format(np.sum(np.abs(F))))
                print('NR number of iterations: {}'.format(itercnt))
                return V, converged

            """ Calculate next voltage guess using NR """

            # matrix pruning and calc of its density factor
            # TODO: need to reimplement the pruning algorithm so that it doesnt
            # create cases of singular jacobians. currently the density factor
            # will mostly choose the non sparse method
            # J[J < 1e-15] = 0.
            # for i in range(0, self.W):
            #     for j in range(0, self.W):
            #         if i != j and J[i,j] < 1e-15:
            #             J[i,j] = 0.
            density = np.count_nonzero(J==0) / np.prod(J.shape)

            # currently just using standard LU factorization
            if density > 1:
                Jn = scipy.sparse.csc_matrix(J)
                lu = scipy.sparse.linalg.splu(Jn)
                dVnew = lu.solve(F)
            else:
                lu, piv = scipy.linalg.lu_factor(J)
                dVnew = scipy.linalg.lu_solve((lu, piv), F)

            # any element in dV is nan, do not update reuse last decrement
            if any(np.isnan(x) for x in dVnew.flatten()):
               pass
            else:
                dV = dVnew

            V = V - dV

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
                if dev.itype == 'ac':
                    iac = dev.ac * np.exp(1j * dev.phase)

                    if self.num_tones == 1 and dev.freq == self.f1:
                        fidx = 1
                    elif self.num_tones == 2:
                        if dev.freq == self.f1:
                            fidx = 4 * self.K2 + 1
                        elif dev.freq == self.f2:
                            fidx = 1
                    else:
                        print('ERROR: only one and two-tones are supported')

                    n1 = dev.n1 - 1
                    if n1 >= 0:
                        Is[S*n1+fidx] += iac.real
                        Is[S*n1+fidx+1] += iac.imag

                    n2 = dev.n2 - 1
                    if n2 >= 0:
                        Is[S*n2+fidx] -= iac.real
                        Is[S*n2+fidx+1] -= iac.imag

                elif dev.itype == 'dc':
                    n1 = dev.n1 - 1
                    if n1 >= 0:
                        Is[S*n1] += dev.dc

                    n2 = dev.n2 - 1
                    if n2 >= 0:
                        Is[S*n2] -= dev.dc

        return Is

    def calc_Y(self):
        # add AC stamps for linear devices
        Y = np.zeros((self.W, self.W))
        for dev in self.lin_devs:
            for k in range(0, self.K+1):
                f = np.abs(self.freqs[k])
                if hasattr(dev, 'add_mthb_stamps') and callable(dev.add_mthb_stamps):
                    dev.add_mthb_stamps(Y, self.S, f, k)

        return Y

    def calc_V0(self):
        # TODO: create better initial solution (maybe TAHB)
        X = DC('HB.DC').run(self.netlist)
        V = np.zeros((self.W,1))
        for n in range(self.N):
            V[n*self.S] = X[n]

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

            # TODO: check for correctness of this
            if self.num_tones == 2:
                for k in range(1,self.K+1):
                    if self.freqs[k] < 0:
                        X[i+2*k] = -X[i+2*k]
        
        return X


