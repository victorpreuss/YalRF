import numpy as np

class Capacitor():
    
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.C    = float(value)

    def get_vsource(self):
        return 0

    def is_linear(self):
        return True

    def add_dc_stamps(self, A, z, x, iidx):
        pass

    def add_ac_stamps(self, A, z, x, iidx, freq):
        y = 1j * 2 * np.pi * freq * self.C
        A[self.n1][self.n1] = A[self.n1][self.n1] + y
        A[self.n2][self.n2] = A[self.n2][self.n2] + y
        A[self.n1][self.n2] = A[self.n1][self.n2] - y
        A[self.n2][self.n1] = A[self.n2][self.n1] - y

    def __str__(self):
        return 'Capacitor: {}\nNodes = {} -> {}\nValue = {}\n'.format(self.name, self.n1, self.n2, self.C)
