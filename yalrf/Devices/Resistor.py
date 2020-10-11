import numpy as np

class Resistor():
    
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.R    = float(value)

    def get_idc(self, x):
        V1 = x[self.n1-1] if self.n1 > 0 else 0.
        V2 = x[self.n2-1] if self.n2 > 0 else 0.
        V = V1 - V2
        return V / self.R

    def get_itran(self, x):
        V = self.get_tran_voltage(x)[:,0]
        return V / self.R

    def get_num_vsources(self):
        return 0

    def is_nonlinear(self):
        return False

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        g = 1.0 / self.R
        A[self.n1][self.n1] = A[self.n1][self.n1] + g
        A[self.n2][self.n2] = A[self.n2][self.n2] + g
        A[self.n1][self.n2] = A[self.n1][self.n2] - g
        A[self.n2][self.n1] = A[self.n2][self.n1] - g

    def add_ac_stamps(self, A, z, x, iidx, freq):
        self.add_dc_stamps(A, z, x, iidx)

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        self.add_dc_stamps(A, z, x, iidx)

    def get_tran_voltage(self, x):
        V1 = x[:,self.n1-1] if self.n1 > 0 else np.zeros((len(x),1))
        V2 = x[:,self.n2-1] if self.n2 > 0 else np.zeros((len(x),1))
        return (V1 - V2)

    def __str__(self):
        return 'Resistor: {}\nNodes = {} -> {}\nValue = {}\n'.format(self.name, self.n1, self.n2, self.R)
