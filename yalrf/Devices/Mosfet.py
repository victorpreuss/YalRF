from scipy.constants import k, e, epsilon_0
import numpy as np

# Mosfet Level 3 model
options = {}
options['Is']   = 1e-14    # bulk junction saturation current
options['N']    = 1.       # bulk junction emission coefficient
options['Vt0']  = 0.7      # zero-bias threshold voltage
options['Lambda'] = 0.     # channel length modulation parameter
options['Kp']   = 2e-5     # transconductance coefficient
options['Gamma'] = 0.      # bulk threshold
options['Phi']  = 0.6      # surface potential
options['Rd']   = 0.       # drain ohmic resistance
options['Rs']   = 0.       # source ohmic resistance
options['Rg']   = 0.       # gate ohmic resistance
options['L']    = 100e-6   # channel length
options['Ld']   = 0.       # lateral diffusion length
options['W']    = 100e-6   # channel width
options['Tox']  = 0.1e-6   # oxide thickness
options['Cgso'] = 0.       # gate-source overlap capacitance per meter of channel width
options['Cgdo'] = 0.       # gate-drain overlap capacitance per meter of channel width
options['Cgbo'] = 0.       # gate-bulk overlap capacitance per meter of channel width
options['Cbd']  = 0.       # zero-bias bulk-drain junction capactiance
options['Cbs']  = 0.       # zero-bias bulk-source junction capactiance
options['Pb']   = 0.8      # bulk junction potential
options['Mj']   = 0.5      # bulk junction bottom grading coefficient
options['Fc']   = 0.5      # bulk junction forward-bias depletion capacitance coefficient
options['Cjsw'] = 0.       # zero-bias junction periphery capacitance per meter of junction perimeter
options['Mjsw'] = 0.33     # bulk junction periphery grading coefficient
options['Tf']   = 0.       # bulk transit time
options['Kf']   = 0.       # flicker noise coefficient
options['Af']   = 1.       # flicker noise exponent
options['Ffe']  = 1.       # flicker noise frequency exponent
options['Nsub'] = 0.       # substrate/bulk doping density
options['Nss']  = 0.       # surface state density
options['Tpg']  = 1.       # gate material type (0=alumina, -1=same as bulk, 1=opposite to bulk)
options['Uo']   = 600.     # surface mobility
options['Rsh']  = 0.       # drain and source diffusion sheet resistance
options['Nrd']  = 1.       # number of equivalent drain squares
options['Nrs']  = 1.       # number of equivalent source squares
options['Cj']   = 0.       # zero-bias bulk junction bottom capacitance per square meter of junction area
options['Js']   = 0.       # bulk junction saturation current per square meter of junction area
options['Ad']   = 0.       # drain diffusion area
options['As']   = 0.       # source diffusion area
options['Pd']   = 0.       # drain junction perimeter
options['Ps']   = 0.       # source junction perimeter
options['Temp'] = 300.     # device temperature
options['Tnom'] = 300.     # parameter measurement temperature

class Mosfet():

    def __init__(self, name, n1, n2, n3, n4=0, ispch=False):
        self.name = name
        self.n1   = n1 # gate (G)
        self.n2   = n2 # drain (D)
        self.n3   = n3 # source (S)
        self.n4   = n4 # bulk (B)

        self.type = -1 if ispch else 1

        self.option = options.copy()
        self.oppoint = {}
        self.adjusted_options = {}

    def get_num_vsource(self):
        return 0

    def is_nonlinear(self):
        return True

    def init(self):
        self.oppoint = {}

        A = self.options['Area']

    def add_dc_stamps(self, A, z, x, iidx):
        self.calc_dc(x)

        G = self.n1
        D = self.n2
        S = self.n3
        B = self.n4
        Vgs = self.oppoint['Vgs']
        Vds = self.oppoint['Vds']
        Vbs = self.oppoint['Vbs']
        Vbd = self.oppoint['Vbd']
        Id  = self.oppoint['Id'] 
        gds = self.oppoint['gds']
        gm  = self.oppoint['gm'] 
        gmb = self.oppoint['gmb']
        Ibd = self.oppoint['Ibd']
        gbd = self.oppoint['gbd']
        Ibs = self.oppoint['Ibs']
        gbs = self.oppoint['gbs']

        Ibdeq = Ibd - gbd * Vbd
        Ibseq = Ibs - gbs * Vbs
        Idseq = Id - gm * Vgs - gmb * Vbs - gds * Vds

        A[G][G] = A[G][G] + 0
        A[G][D] = A[G][D] + 0
        A[G][S] = A[G][S] + 0
        A[G][B] = A[G][B] + 0
        A[D][G] = A[D][G] + gm
        A[D][D] = A[D][D] + gds + gbd
        A[D][S] = A[D][S] - gds - gm - gmb
        A[D][B] = A[D][B] + gmb - gbd
        A[S][G] = A[S][G] - gm
        A[S][D] = A[S][D] - gds
        A[S][S] = A[S][S] + gbs + gds + gm + gmb
        A[S][B] = A[S][B] - gbs - gmb
        A[B][G] = A[B][G] + 0
        A[B][D] = A[B][D] - gbd
        A[B][S] = A[B][S] - gbs
        A[B][B] = A[B][B] + gbs + gbd
        z[G] = z[G] + 0
        z[D] = z[D] + (+ Ibdeq - Idseq) * self.type
        z[S] = z[S] + (+ Ibseq + Idseq) * self.type
        z[B] = z[B] + (- Ibdeq - Ibseq) * self.type

    def add_ac_stamps(self, A, z, x, iidx, freq):
        pass

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        pass

    def calc_oppoint(self, x, usevlimit=True):
        pass

    def calc_dc(self, x, usevlimit=True):
        G   = self.n1
        D   = self.n2
        S   = self.n3
        B   = self.n4

        Vt  = k * self.options['Temp'] / q
        L   = self.options['L']
        Ld  = self.options['Ld']
        W   = self.options['W']
        Kp  = self.options['Kp']
        Vt0 = self.options['Vt0']
        Phi = self.options['Phi']
        Gamma = self.options['Gamma']
        Lambda = self.options['Lambda']
        Is  = self.options['Is']
        N   = self.options['N']

        Vg = x[G-1,0] if G > 0 else 0
        Vd = x[D-1,0] if D > 0 else 0
        Vs = x[S-1,0] if S > 0 else 0
        Vb = x[B-1,0] if B > 0 else 0

        Vgs = (Vg - Vs) * self.type
        Vds = (Vd - Vs) * self.type
        Vbs = (Vb - Vs) * self.type
        Vbd = (Vb - Vd) * self.type

        # TODO: include here limiting algorithms for the fet and
        #       for the parasitic junction diodes to avoid overflow

        Leff = L - 2 * Ld

        if Tox > 0:
            Eps0 = epsilon_0 # vacuum permittivity
            EpsR = 3.9       # SiO2 relative permittivity
            Cox = EpsR * Eps0 / Tox
        else:
            Cox = 0

        if Kp > 0:
            Beta = Kp * W / Leff
        else:
            Beta = 1e-4 * Uo * Cox * W / Leff

        if Vds >= 0:
            Vth = Vt0 + Gamma * (np.sqrt(Phi-Vbs) - np.sqrt(Phi))
        else:
            Vth = Vt0 + Gamma * (np.sqrt(Phi-Vbd) - np.sqrt(Phi))

        if Vgs - Vth <= 0:
            Id  = 0
            gds = 0
            gm  = 0
            gmb = 0
        elif Vgs - Vth > 0 and Vgs - Vth < Vds:
            Id  = Beta / 2 * (1 + Lambda * Vds) * np.square(Vgs - Vth)
            gds = Beta / 2 * Lambda * np.square(Vgs - Vth)
            gm  = Beta * (1 + Lambda * Vds) * np.square(Vgs - Vth)
            gmb = gm * Gamma / (2 * np.sqrt(Phi - Vbs))
        elif Vgs - Vth > Vds:
            Id  = Beta * (1 + Lambda * Vds) * (Vgs - Vth - Vds / 2 ) * Vds
            gds = Beta * (1 + Lambda * Vds) * (Vgs - Vth - Vds) + Beta * Lambda * Vds * (Vgs - Vth - Vds / 2)
            gm  = Beta * (1 + Lambda * Vds) * Vds
            gmb = gm * Gamma / (2 * np.sqrt(Phi - Vbs))

        Ibd = Is * np.exp(Vbd / (N * Vt) - 1)
        gbd = Is / (N * Vt) * np.exp(Vbd / (N * Vt) - 1)

        Ibs = Is * np.exp(Vbs / (N * Vt) - 1)
        gbs = Is / (N * Vt) * np.exp(Vbs / (N * Vt) - 1)

        self.oppoint['Vgs'] = Vgs
        self.oppoint['Vds'] = Vds
        self.oppoint['Vbs'] = Vbs
        self.oppoint['Vbd'] = Vbd
        self.oppoint['Id']  = Id
        self.oppoint['gds'] = gds
        self.oppoint['gm']  = gm
        self.oppoint['gmb'] = gmb
        self.oppoint['Ibd'] = Ibd
        self.oppoint['gbd'] = gbd
        self.oppoint['Ibs'] = Ibs
        self.oppoint['gbs'] = gbs


