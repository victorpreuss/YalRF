import numpy as np

class Inductor():
    
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.L    = float(value)

    def get_num_vsources(self, analysis):
        if analysis == 'dc':
            return 1
        elif analysis == 'ac':
            return 0

    def is_linear(self):
        return True

    def add_dc_stamps(self, A, z, x, iidx):
        A[self.n1][iidx] = +1.0 
        A[self.n2][iidx] = -1.0 
        A[iidx][self.n1] = +1.0 
        A[iidx][self.n2] = -1.0
        z[iidx] = 0.0

    def add_ac_stamps(self, A, z, x, iidx, freq):
        y = 1.0 / (1j * 2 * np.pi * freq * self.L)
        A[self.n1][self.n1] = A[self.n1][self.n1] + y
        A[self.n2][self.n2] = A[self.n2][self.n2] + y
        A[self.n1][self.n2] = A[self.n1][self.n2] - y
        A[self.n2][self.n1] = A[self.n2][self.n1] - y

    def __str__(self):
        return 'Inductor: {}\nNodes = {} -> {}\nValue = {}\n'.format(self.name, self.n1, self.n2, self.L)
