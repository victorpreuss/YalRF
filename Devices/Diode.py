from scipy.constants import k, e
import numpy as np

# diode model options
options = dict()
options['Temp'] = 300.0 # model temperature
options['Is'] = 1e-15   # saturation current
options['N']  = 1.0     # emission coefficient

class Diode():
    
    def __init__(self, name, n1, n2):
        self.name = name
        self.n1   = n1 # anode (+)
        self.n2   = n2 # cathode (-)

        self.options = options.copy() # diode options

        # useful for limiting scheme
        self.Vdold = 0

    def get_num_vsources(self, analysis):
        return 0

    def is_linear(self):
        return False

    def add_dc_stamps(self, A, z, x, iidx):

        Vt = k * self.options['Temp'] / e
        Is = self.options['Is']
        N  = self.options['N']

        V1 = x[self.n1-1] if self.n1 else 0.0
        V2 = x[self.n2-1] if self.n2 else 0.0

        Vd = V1 - V2

        # limiting algorithm to avoid overflow
        Vcrit = N * Vt * np.log(N * Vt / (np.sqrt(2) * Is))
        if (Vd > Vcrit and Vd > 0):
            Vd = self.Vdold + N * Vt * np.log(1 + (Vd - self.Vdold) / (N * Vt))
        self.Vdold = Vd

        Id = Is * (np.exp(Vd / (N * Vt)) - 1)

        g  = Is / (N * Vt) * np.exp(Vd / (N * Vt))
        I  = Id - g * Vd

        A[self.n1][self.n1] = A[self.n1][self.n1] + g
        A[self.n2][self.n2] = A[self.n2][self.n2] + g
        A[self.n1][self.n2] = A[self.n1][self.n2] - g
        A[self.n2][self.n1] = A[self.n2][self.n1] - g
        z[self.n1] = z[self.n1] - I
        z[self.n2] = z[self.n2] + I

    def __str__(self):
        return 'Diode: {}\nNodes = {} -> {}\n'.format(self.name, self.n1, self.n2)
