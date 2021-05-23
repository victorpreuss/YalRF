import numpy as np

class Inductor():
    
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.L    = float(value)

    def get_num_vsources(self):
        return 1

    def is_nonlinear(self):
        return False

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        A[self.n1][iidx] = +1.0
        A[self.n2][iidx] = -1.0
        A[iidx][self.n1] = +1.0
        A[iidx][self.n2] = -1.0
        z[iidx] = 0.0

    def add_ac_stamps(self, A, z, x, iidx, freq):
        y = 1.0 / (1j * 2 * np.pi * freq * self.L) # if freq != 0. else 1e9
        # A[self.n1][self.n1] += y
        # A[self.n2][self.n2] += y
        # A[self.n1][self.n2] -= y
        # A[self.n2][self.n1] -= y
        A[self.n1][iidx] = +1.0
        A[self.n2][iidx] = -1.0
        A[iidx][self.n1] = +y
        A[iidx][self.n2] = -y
        A[iidx][iidx] = -1.0

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        # get inductor voltage and current
        V1 = xt[-1][self.n1-1] if self.n1 > 0 else 0.
        V2 = xt[-1][self.n2-1] if self.n2 > 0 else 0.
        Vn = V1 - V2
        In = xt[-1][iidx-1]

        # calculate companion model parameters

        # trapezoidal
        geq = tstep / (2. * self.L)
        Ieq = In + geq * Vn

        # implicit euler
        # geq = tstep / self.L
        # Ieq = In

        # add to MNA
        A[self.n1][iidx] = +1.0
        A[self.n2][iidx] = -1.0
        A[iidx][self.n1] = +geq
        A[iidx][self.n2] = -geq
        A[iidx][iidx] = -1.0
        z[iidx] = z[iidx] - Ieq

    def add_mthb_stamps(self, Y, S, freq, freqidx):
        if freqidx == 0:
            n = (self.n1 - 1) * S
            m = (self.n2 - 1) * S
            y = 1e6
            if (self.n1 > 0):
                Y[n,n] += y
            if (self.n2 > 0):
                Y[m,m] += y
            if (self.n1 > 0 and self.n2 > 0):
                Y[n,m] -= y
                Y[m,n] -= y
        else:
            n = (self.n1 - 1) * S + 2 * (freqidx - 1) + 1
            m = (self.n2 - 1) * S + 2 * (freqidx - 1) + 1
            y = 1. / (1j * 2 * np.pi * freq * self.L)
            Ymnk = np.array([[y.real, -y.imag],[y.imag, y.real]])
            if (self.n1 > 0):
                Y[n:n+2,n:n+2] += Ymnk 
            if (self.n2 > 0):
                Y[m:m+2,m:m+2] += Ymnk 
            if (self.n1 > 0 and self.n2 > 0):
                Y[n:n+2,m:m+2] -= Ymnk 
                Y[m:m+2,n:n+2] -= Ymnk 

    def __str__(self):
        return 'Inductor: {}\nNodes = {} -> {}\nValue = {}\n'.format(self.name, self.n1, self.n2, self.L)
