from scipy.constants import k, e
import numpy as np

# diode model options
options = {}
options['Temp'] = 300.0  # device temperature
options['Is']  = 1e-15   # saturation current
options['N']   = 1.0     # emission coefficient
options['Isr'] = 0.0     # recombination current parameter
options['Nr']  = 2.0     # emission coefficient for Isr
options['Rs']  = 0.0     # ohmic resistance
options['Cj0'] = 0.0     # zero-bias junction capacitance
options['M']   = 0.5     # grading coefficient
options['Vj']  = 0.7     # junction potential
options['Fc']  = 0.5     # forward-bias depletion capacitance coefficient
options['Cp']  = 0.0     # linear capacitance
options['Tt']  = 0.0     # transit time
options['Bv']  = 1e15    # reverse breakdown voltage
options['Ibv'] = 0.001   # current at reverse breakdown voltage
options['Ikf'] = 1e12    # high-injection knee current
options['Kf']  = 0.0     # flicker noise coefficient
options['Af']  = 1.0     # flicker noise exponent
options['Ffe'] = 1.0     # flicker noise frequency exponent
options['Xti'] = 3.0     # saturation current exponent
options['Eg']  = 1.11    # energy bandgap
options['Ega'] = 7.02e-4 # temperature coefficient 1 for energy bandgap
options['Egb'] = 1108.0  # temperature coefficient 2 for energy bandgap
options['Tbv'] = 0.0     # Bv linear temperature coefficient
options['Trs'] = 0.0     # Rs linear temperature coefficient
options['Tt1'] = 0.0     # Tt linear temperature coefficient
options['Tt2'] = 0.0     # Tt quadratic temperature coefficient
options['Tm1'] = 0.0     # M linear temperature coefficient
options['Tm2'] = 0.0     # M quadratic temperature coefficient
options['Tnom'] = 300.0  # temperature at which parameters were extracted
options['Area'] = 1.0    # diode area multiplier

# TODO: implement reverse breakdown, temperature dependence and area dependence
class Diode():
    
    def __init__(self, name, n1, n2):
        self.name = name
        self.n1   = n1 # anode (+)
        self.n2   = n2 # cathode (-)

        self.options = options.copy() # diode options

        # for the limiting scheme
        self.Vdold = 0.

    def get_num_vsources(self, analysis):
        return 0

    def is_linear(self):
        return False

    def add_dc_stamps(self, A, z, x, iidx):
        T = self.options['Temp']
        Vt = k * T / e

        V1 = x[self.n1-1] if self.n1 != 0 else 0.
        V2 = x[self.n2-1] if self.n2 != 0 else 0.
        Vd = V1 - V2

        Vd = self.limit_diode_voltage(Vd, Vt)
        gd, Id = self.get_gd_and_Id(Vd, Vt)
        I = Id - gd * Vd

        A[self.n1][self.n1] = A[self.n1][self.n1] + gd
        A[self.n2][self.n2] = A[self.n2][self.n2] + gd
        A[self.n1][self.n2] = A[self.n1][self.n2] - gd
        A[self.n2][self.n1] = A[self.n2][self.n1] - gd
        z[self.n1] = z[self.n1] - I
        z[self.n2] = z[self.n2] + I

    def add_ac_stamps(self, A, z, x, iidx, freq):
        Cj0 = self.options['Cj0']
        M = self.options['M']
        Vj = self.options['Vj']
        Fc = self.options['Fc']
        Cp = self.options['Cp']
        Tt = self.options['Tt']
        T = self.options['Temp']
        Vt = k * T / e

        V1 = x[self.n1-1] if self.n1 != 0 else 0.
        V2 = x[self.n2-1] if self.n2 != 0 else 0.
        Vd = V1 - V2

        Vd = self.limit_diode_voltage(Vd, Vt)
        gd, Id = self.get_gd_and_Id(Vd, Vt)

        if Vd / Vj <= Fc:
            Cj = Cj0 * np.power((1. - Vd / Vj), -M)
            Qj = Cj0 * Vj / (1. - M) * (1. - np.power((1. - Vd / Vj), (1. - M)))
        else:
            Cj = Cj0 / np.power((1. - Fc), M) * (1. + M * (Vd / Vj - Fc) / (1. - Fc))
            Xd = (1. - np.power((1. - Fc), (1. - M))) / (1. - M) + \
                 (1. - Fc * (1. + M)) / np.power((1. - Fc), (1. + M)) * (Vd / Vj - Fc) + \
                 M / (2. * np.power((1. - Fc), (1. + M))) * (np.square(Vd / Vj) - np.square(Fc))
            Qj = Cj0 * Vj * Xd
        
        Cd = Cp + Tt * gd + Cj
        Qd = Cp * Vd + Tt * Id + Qj
        y = gd + 1j * 2. * np.pi * freq * Cd
        
        # print('np.pi = {}\nfreq = {}\n'.format(np.pi, freq))
        # print('Vd = {}\nId = {}\ngd = {}\nCd = {}\ny = {}\n'.format(Vd, Id, gd, Cd, y))

        A[self.n1][self.n1] = A[self.n1][self.n1] + y
        A[self.n2][self.n2] = A[self.n2][self.n2] + y
        A[self.n1][self.n2] = A[self.n1][self.n2] - y
        A[self.n2][self.n1] = A[self.n2][self.n1] - y

    def limit_diode_voltage(self, Vd, Vt):
        Is = self.options['Is']
        N = self.options['N']
        Vcrit = N * Vt * np.log(N * Vt / (np.sqrt(2.) * Is))
        if (Vd > 0. and Vd > Vcrit):
            Vd = self.Vdold + N * Vt * np.log1p((Vd - self.Vdold) / (N * Vt))
        self.Vdold = Vd
        return Vd

    def get_gd_and_Id(self, Vd, Vt):
        Is = self.options['Is']
        N = self.options['N']
        Isr = self.options['Isr']
        Nr = self.options['Nr']
        Ikf = self.options['Ikf']

        # forward current
        Idf = Is * np.expm1(Vd / (N * Vt))
        gdf = Is / (N * Vt) * np.exp(Vd / (N * Vt))

        # reverse current
        Idr = Isr * np.expm1(Vd / (Nr * Vt))
        gdr = Isr / (Nr * Vt) * np.exp(Vd / (Nr * Vt))

        # high injection
        if Ikf > 0.:
            Idf = Idf * np.sqrt(Ikf / (Ikf + Idf))
            gdf = gdf * (1. - 0.5 * Idf / (Ikf + Idf)) * np.sqrt(Ikf / (Ikf + Idf))

        # diode total current
        Id = Idf + Idr
        gd = gdf + gdr

        return gd, Id

    # a simple exponential diode, used for experiments :-)
    def get_gd_and_Id_shockley(self, A, z, x, iidx):
        Is = self.options['Is']
        N  = self.options['N']
        Vt = k * self.options['Temp'] / e

        V1 = x[self.n1-1] if self.n1 else 0.0
        V2 = x[self.n2-1] if self.n2 else 0.0

        Vd = V1 - V2

        # limiting algorithm to avoid overflow
        Vcrit = N * Vt * np.log(N * Vt / (np.sqrt(2.) * Is))
        if (Vd > Vcrit and Vd > 0.):
            Vd = self.Vdold + N * Vt * np.log(1. + (Vd - self.Vdold) / (N * Vt))
        self.Vdold = Vd

        Id = Is * (np.exp(Vd / (N * Vt)) - 1.)

        g = Is / (N * Vt) * np.exp(Vd / (N * Vt))
        I = Id - g * Vd

        return g, I

    def __str__(self):
        return 'Diode: {}\nNodes = {} -> {}\n'.format(self.name, self.n1, self.n2)
