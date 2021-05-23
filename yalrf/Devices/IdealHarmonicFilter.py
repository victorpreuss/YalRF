import numpy as np

class IdealHarmonicFilter():
    
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.freq = float(value)

        # conductance of 'shorted' harmonics
        self.g = 1e9

    def get_num_vsources(self):
        return 0

    def is_nonlinear(self):
        return False

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        pass

    def add_ac_stamps(self, A, z, x, iidx, freq):
        g = self.g if freq == self.freq else 0.
        A[self.n1][self.n1] = A[self.n1][self.n1] + g
        A[self.n2][self.n2] = A[self.n2][self.n2] + g
        A[self.n1][self.n2] = A[self.n1][self.n2] - g
        A[self.n2][self.n1] = A[self.n2][self.n1] - g

    def add_mthb_stamps(self, Y, S, freq, freqidx):
        gmin = 1e-12 # conductance of open harmonics
        if freqidx == 0:
            n = (self.n1 - 1) * S
            m = (self.n2 - 1) * S
            y = gmin
            if (self.n1 > 0):
                Y[n,n] += y
            if (self.n2 > 0):
                Y[m,m] += y
            if (self.n1 > 0 and self.n2 > 0):
                Y[n,m] -= y
                Y[m,n] -= y
        else:
            y = self.g if freq == self.freq else gmin
            n = (self.n1 - 1) * S + 2 * (freqidx - 1) + 1
            m = (self.n2 - 1) * S + 2 * (freqidx - 1) + 1
            Ymnk = np.array([[y, 0],[0, y]])
            if (self.n1 > 0):
                Y[n:n+2,n:n+2] += Ymnk 
            if (self.n2 > 0):
                Y[m:m+2,m:m+2] += Ymnk 
            if (self.n1 > 0 and self.n2 > 0):
                Y[n:n+2,m:m+2] -= Ymnk 
                Y[m:m+2,n:n+2] -= Ymnk 

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        self.add_dc_stamps(A, z, x, iidx)

    def __str__(self):
        return 'IdealHarmonicFilter: {}\nNodes = {} -> {}\nValue = {}\n'.format(self.name, self.n1, self.n2, self.freq)
