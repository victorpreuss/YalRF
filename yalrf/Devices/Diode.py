from math import inf, isinf, isfinite
from scipy.constants import k, e
import numpy as np

# Diode model options
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
options['Bv']  = inf     # reverse breakdown voltage
options['Ibv'] = 0.001   # current at reverse breakdown voltage
options['Ikf'] = inf     # high-injection knee current
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

# TODO: temperature dependence
#       noise
class Diode():
    
    def __init__(self, name, n1, n2):
        self.name = name
        self.n1   = n1 # anode (+)
        self.n2   = n2 # cathode (-)

        self.options = options.copy() # diode options
        self.oppoint = {}

        # this options vector holds the corrected values
        # for the diode according to area and temperature
        self.adjusted_options = {}

        # for the limiting scheme
        self.Vdold = 0.

    def get_num_vsources(self, analysis):
        return 0

    def is_nonlinear(self):
        return True

    def init(self):
        # clear model internal variables
        self.oppoint = {}
        self.Vdold = 0.
        
        # area and temperature dependent adjustments
        A = self.options['Area']
        self.adjusted_options['Is'] = self.options['Is'] * A
        self.adjusted_options['Isr'] = self.options['Isr'] * A
        self.adjusted_options['Ikf'] = self.options['Ikf'] * A
        self.adjusted_options['Ibv'] = self.options['Ibv'] * A
        self.adjusted_options['Cj0'] = self.options['Cj0'] * A
        self.adjusted_options['Rs'] = self.options['Rs'] / A

    def add_dc_stamps(self, A, z, x, iidx):
        # calculate dc parameters
        self.calc_dc(x)

        It = self.oppoint['It']
        gt = self.oppoint['gt']

        A[self.n1][self.n1] = A[self.n1][self.n1] + gt
        A[self.n2][self.n2] = A[self.n2][self.n2] + gt
        A[self.n1][self.n2] = A[self.n1][self.n2] - gt
        A[self.n2][self.n1] = A[self.n2][self.n1] - gt
        z[self.n1] = z[self.n1] - It
        z[self.n2] = z[self.n2] + It

    def add_ac_stamps(self, A, z, x, iidx, freq):
        gd = self.oppoint['gd']
        Cd = self.oppoint['Cd']
        Rs = self.adjusted_options['Rs']
        y = 1. / (Rs + 1. / (gd + 1j * 2. * np.pi * freq * Cd))

        A[self.n1][self.n1] = A[self.n1][self.n1] + y
        A[self.n2][self.n2] = A[self.n2][self.n2] + y
        A[self.n1][self.n2] = A[self.n1][self.n2] - y
        A[self.n2][self.n1] = A[self.n2][self.n1] - y

    def calc_oppoint(self, x):
        self.calc_dc(x)

        Cj0 = self.adjusted_options['Cj0']
        M = self.options['M']
        Vj = self.options['Vj']
        Fc = self.options['Fc']
        Cp = self.options['Cp']
        Tt = self.options['Tt']
        Vd = self.oppoint['Vd']
        gd = self.oppoint['gd']

        if Vd / Vj <= Fc:
            Cj = Cj0 * np.power((1. - Vd / Vj), -M)
        else:
            Cj = Cj0 / np.power((1. - Fc), M) * (1. + M * (Vd / Vj - Fc) / (1. - Fc))

        Cd = Cp + Tt * gd + Cj

        self.oppoint['Cd'] = Cd
        self.oppoint['Cj'] = Cj
        self.oppoint['Cdiff'] = Tt * gd

        # if Vd / Vj <= Fc:
        #     Qj = Cj0 * Vj / (1. - M) * (1. - np.power((1. - Vd / Vj), (1. - M)))
        # else:
        #     Xd = (1. - np.power((1. - Fc), (1. - M))) / (1. - M) + \
        #          (1. - Fc * (1. + M)) / np.power((1. - Fc), (1. + M)) * (Vd / Vj - Fc) + \
        #          M / (2. * np.power((1. - Fc), (1. + M))) * (np.square(Vd / Vj) - np.square(Fc))
        #     Qj = Cj0 * Vj * Xd
        # Qd = Cp * Vd + Tt * Id + Qj

    def calc_dc(self, x):
        Is = self.adjusted_options['Is']
        N = self.options['N']
        Isr = self.adjusted_options['Isr']
        Nr = self.options['Nr']
        Ikf = self.adjusted_options['Ikf']
        Bv = self.options['Bv']
        Ibv = self.adjusted_options['Ibv']
        Rs = self.adjusted_options['Rs']
        T = self.options['Temp']
        Vt = k * T / e

        # get diode voltage
        V1 = x[self.n1-1,0] if self.n1 > 0 else 0.
        V2 = x[self.n2-1,0] if self.n2 > 0 else 0.

        # remove voltage drop caused by series resistance
        Vtot = V1 - V2
        if 'It' in self.oppoint:
            It = self.oppoint['It']
            gt = self.oppoint['gt']
            Id = It + gt * Vtot
            Vd = Vtot - Id * Rs
        else:
            Vd = Vtot # ignore Rs effect at first iteration

        Vd = self.limit_diode_voltage(Vd, Vt)

        # forward current
        Idf = Is * (exp_lim(Vd / (N * Vt)) - 1.)
        gdf = Is / (N * Vt) * exp_lim(Vd / (N * Vt))

        # reverse current
        Idr = Isr * (exp_lim(Vd / (Nr * Vt)) - 1.)
        gdr = Isr / (Nr * Vt) * exp_lim(Vd / (Nr * Vt))

        # high injection
        if isfinite(Ikf):
            Idf = Idf * np.sqrt(Ikf / (Ikf + Idf))
            gdf = gdf * (1. - 0.5 * Idf / (Ikf + Idf)) * np.sqrt(Ikf / (Ikf + Idf))

        # reverse breakdown
        Ibr = 0.
        gbr = 0.
        if isfinite(Bv):
            Ibr = Ibv * exp_lim(- Bv / Vt) * (1. - exp_lim(- Vd / Vt))
            gbr = Ibv / Vt * exp_lim(- (Bv + Vd) / Vt)

        # diode total current
        Id = Idf + Idr + Ibr
        gd = gdf + gdr + gbr

        # add effect of series resistance
        It = (Id - gd * Vd) / (1. + gd * Rs)
        gt = gd / (1. + gd * Rs)

        self.oppoint['Vd'] = Vd
        self.oppoint['Id'] = Id
        self.oppoint['gd'] = gd
        self.oppoint['It'] = It
        self.oppoint['gt'] = gt
        self.oppoint['Idf'] = Idf
        self.oppoint['gdf'] = gdf
        self.oppoint['Idr'] = Idr
        self.oppoint['gdr'] = gdr
        self.oppoint['Ibr'] = Ibr
        self.oppoint['gbr'] = gbr

    def limit_diode_voltage(self, Vd, Vt):
        Is = self.adjusted_options['Is']
        N = self.options['N']

        # limiting schemes
        # Vcrit = N * Vt * np.log(N * Vt / (np.sqrt(2.) * Is))
        # if (Vd > 0. and Vd > Vcrit):
        #     Vd = self.Vdold +  N * Vt * np.log1p((Vd - self.Vdold) / (N * Vt))

        # if Vd > 10 * Vt:
        #     if Vd > self.Vdold:
        #         Vd = self.Vdold + 2 * Vt
        #     else:
        #         Vd = self.Vdold - 2 * Vt
        
        Vd = self.Vdold + 10. * N * Vt * np.tanh((Vd - self.Vdold) / (10. * N * Vt))
        self.Vdold = Vd

        return Vd

    # Legacy: a simple exponential diode, used for initial experiments :-)
    def get_gd_and_Id_shockley(self, A, z, x, iidx):
        Is = self.adjusted_options['Is']
        N  = self.options['N']
        Vt = k * self.options['Temp'] / e

        V1 = x[self.n1-1] if self.n1 else 0.
        V2 = x[self.n2-1] if self.n2 else 0.

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

# limit the maximum derivative of the exponential function
def exp_lim(x):
    return np.exp(x) if x < 70. else np.exp(70.) + np.exp(70.) * (x - 70.)
