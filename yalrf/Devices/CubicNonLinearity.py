import numpy as np

class CubicNonLinearity():
    
    def __init__(self, name, n1, n2, alpha=1e-3):
        self.name  = name
        self.n1    = n1
        self.n2    = n2
        self.alpha = alpha

    def get_num_vsources(self):
        return 0

    def is_nonlinear(self):
        return True

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        V1 = x[self.n1-1] if self.n1 > 0 else 0.
        V2 = x[self.n2-1] if self.n2 > 0 else 0.
        V = V1 - V2

        g = 2 * self.alpha * V*V
        I = self.alpha * V*V*V
        Ieq = I - g * V

        A[self.n1][self.n1] = A[self.n1][self.n1] + g
        A[self.n2][self.n2] = A[self.n2][self.n2] + g
        A[self.n1][self.n2] = A[self.n1][self.n2] - g
        A[self.n2][self.n1] = A[self.n2][self.n1] - g
        z[self.n1] = z[self.n1] - Ieq
        z[self.n2] = z[self.n2] + Ieq

    def get_mthb_params(self, V, Vold):
        g = 2 * self.alpha * V*V
        I = self.alpha * V*V*V
        Ieq = I - g * V

        return g, I

    def __str__(self):
        return 'CubicNL: {}\nNodes = {} -> {}\n'.format(self.name, self.n1, self.n2)

