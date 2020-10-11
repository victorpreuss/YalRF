import numpy as np

class Capacitor():
    
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.C    = float(value)

        self.I = [] # list of accepted computed currents during transient

    def get_idc(self, x):
        return 0.

    def get_itran(self, x):
        return np.array(self.I)

    def get_num_vsources(self):
        return 0

    def is_nonlinear(self):
        return False

    def init(self):
        # initial current at capacitor is zero
        self.I = [0.]

    def add_dc_stamps(self, A, z, x, iidx):
        pass

    def add_ac_stamps(self, A, z, x, iidx, freq):
        y = 1j * 2 * np.pi * freq * self.C
        A[self.n1][self.n1] = A[self.n1][self.n1] + y
        A[self.n2][self.n2] = A[self.n2][self.n2] + y
        A[self.n1][self.n2] = A[self.n1][self.n2] - y
        A[self.n2][self.n1] = A[self.n2][self.n1] - y

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        # get capacitor voltage and current
        Vn = self.get_voltage(xt[-1])
        In = self.I[-1]

        # calculate companion model parameters

        # trapezoidal
        geq = 2. * self.C / tstep
        Ieq = - In - geq * Vn

        # implicit euler
        # self.geq = self.C / tstep
        # self.Ieq = - self.geq * Vn

        # add to MNA
        A[self.n1][self.n1] = A[self.n1][self.n1] + geq
        A[self.n2][self.n2] = A[self.n2][self.n2] + geq
        A[self.n1][self.n2] = A[self.n1][self.n2] - geq
        A[self.n2][self.n1] = A[self.n2][self.n1] - geq
        z[self.n1] = z[self.n1] - Ieq
        z[self.n2] = z[self.n2] + Ieq

    def save_tran(self, xt, tstep):
        # get last capacitor voltages
        Vnprev = self.get_voltage(xt[-2])
        Vn = self.get_voltage(xt[-1])

        # calculate capacitor current for storage
        geq = (2. * self.C / tstep)
        In = geq * (Vn - Vnprev) - self.I[-1]
        self.I.append(In)

    def get_voltage(self, x):
        V1 = x[self.n1-1] if self.n1 > 0 else 0.
        V2 = x[self.n2-1] if self.n2 > 0 else 0.
        return (V1 - V2)

    def __str__(self):
        return 'Capacitor: {}\nNodes = {} -> {}\nValue = {}\n'.format(self.name, self.n1, self.n2, self.C)
