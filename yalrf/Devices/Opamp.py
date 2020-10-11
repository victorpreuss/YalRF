import numpy as np

class Opamp():
    
    def __init__(self, name, n1, n2, n3, G=100e3, Vmax=10e3):
        self.name = name
        self.n1   = n1 # (+)
        self.n2   = n2 # (-)
        self.n3   = n3 # out

        self.G    = float(G)
        self.Vmax = float(Vmax)

        self.oppoint = {}

    def get_num_vsources(self):
        return 1

    def is_nonlinear(self):
        return True

    def init(self):
        self.oppoint = {}

    def add_dc_stamps(self, A, z, x, iidx):
        # calculate dc parameters
        self.calc_dc(x)

        g = self.oppoint['g']
        Veq = self.oppoint['Veq']

        A[iidx][self.n1] = +g
        A[iidx][self.n2] = -g
        A[self.n3][iidx] = +1.
        A[iidx][self.n3] = -1.
        z[iidx] = Veq

    def add_ac_stamps(self, A, z, x, iidx, freq):
        g = self.oppoint['g']
        A[iidx][self.n1] = +g
        A[iidx][self.n2] = -g
        A[self.n3][iidx] = +1.
        A[iidx][self.n3] = -1.

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        self.add_dc_stamps(A, z, x, iidx)

    def calc_dc(self, x):
        Vin, Vout = self.get_vdc(x)
        g = self.G / (1. + np.square(np.pi / (2. * self.Vmax) * self.G * Vin))
        Veq = g * Vin - self.Vmax * 2. / np.pi * np.arctan(np.pi / (2. * self.Vmax) * self.G * Vin)

        self.oppoint['g'] = g
        self.oppoint['Veq'] = Veq

    def calc_oppoint(self, x):
        self.calc_dc(x)

    def save_oppoint(self):
        pass

    def save_tran(self, x, tstep):
        pass

    def get_vdc(self, x):
        V1 = x[self.n1-1] if self.n1 > 0 else 0.
        V2 = x[self.n2-1] if self.n2 > 0 else 0.
        Vin = V1 - V2
        Vout = x[self.n3-1] if self.n2 > 0 else 0.
        return Vin, Vout

    def __str__(self):
        return 'Opamp: {}\nNodes = {}, {}, {}\nValue = {}\n'.format(self.name, self.n1, self.n2, self.n3, self.G)
