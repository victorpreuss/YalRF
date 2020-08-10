from scipy.constants import k, e
import numpy as np

# BJT Spice Gummel-Poon model options
options = {}
options['Temp'] = 300.0 # device temperature
options['Is']  = 1e-15  # saturation current
options['Nf']  = 1.0    # forward emission coefficient
options['Nr']  = 1.0    # reverse emission coefficient
options['Ikf'] = 1e12   # high current corner for forward beta
options['Ikr'] = 1e12   # high current corner for reverse beta
options['Vaf'] = 1e15   # forward Early Voltage
options['Var'] = 1e15   # reverse Early Voltage
options['Ise'] = 0.0    # base-emitter leakage saturation current
options['Ne']  = 1.5    # base-emitter leakage emission coefficient
options['Isc'] = 0.0    # base-collector leakage saturation current
options['Nc']  = 2.0    # base-collector leakage emission coefficient
options['Bf']  = 100.0  # forward beta
options['Br']  = 1.0    # reverse beta
options['Rbm'] = 0.0    # minimum base resistance for high-currents
options['Irb'] = 1e12   # current for base resistance midpoint
options['Rc']  = 0.0    # collector ohmic resistance
options['Re']  = 0.0    # emitter ohmic resistance
options['Rb']  = 0.0    # zero-bias base resistance
options['Cje'] = 0.0    # base-emitter zero-bias depletion capacitance
options['Vje'] = 0.75   # base-emitter junction built-in potential
options['Mje'] = 0.33   # base-emitter junction exponential factor
options['Cjc'] = 0.0    # base-collector zero-bias depletion capacitance
options['Vjc'] = 0.75   # base-collector junction built-in potential
options['Mjc'] = 0.33   # base-collector junction exponential factor
options['Xcjc'] = 1.0   # fraction of Cjc that goes to internal base pin
options['Cjs'] = 0.0    # zero-bias collector-substrate capacitance
options['Vjs'] = 0.75   # substrate junction built-in potential
options['Mjs'] = 0.0    # substrate junction exponential factor
options['Fc']  = 0.5    # forward-bias depletion capacitance coefficient
options['Tf']  = 0.0    # ideal forward transit time
options['Xtf'] = 0.0    # coefficient of bias-dependence for Tf
options['Vtf'] = 1e15   # coefficient of base-collector voltage dependence for Tf
options['Itf'] = 0.0    # high-current effect on Tf
options['Ptf'] = 0.0    # excess phase at frequency 1 / (2 pi Tf)
options['Tr']  = 0.0    # ideal reverse transit time
options['Kf']  = 0.0    # flicker noise coefficient
options['Af']  = 1.0    # flicker noise exponent
options['Ffe'] = 1.0    # flicker noise frequency exponent
options['Kb']  = 0.0    # burst noise coefficient
options['Ab']  = 1.0    # burst noise exponent
options['Fb']  = 1.0    # burst noise corner frequency
options['Xti'] = 3.0    # saturation current exponent
options['Xtb'] = 0.0    # temperature exponent for Bf and Br
options['Eg']  = 1.11   # energy bandgap
options['Tnom'] = 300.0 # temperature at which parameters were extracted
options['Area'] = 1.0   # bjt area multiplier

# TODO: base resistance (bias dependent)
#       DC resistances Rb, Rc and Re
#       excess phase
#       temperature dependence
#       area dependence
#       Cbcxdep is currently unused
#       noise
class BJT():
    
    def __init__(self, name, n1, n2, n3, n4=0):
        self.name = name
        self.n1   = n1 # base (B)
        self.n2   = n2 # collector (C)
        self.n3   = n3 # emitter (E)
        self.n4   = n4 # sbustrate (S)

        self.options = options.copy() # bjt options
        self.oppoint = {}

        # this options vector holds the corrected values
        # for the BJT according to area and temperature
        self.adjusted_options = options.copy()

        # for the limiting scheme
        self.Vbeold = 0
        self.Vbcold = 0

    def get_num_vsources(self, analysis):
        return 0

    def is_nonlinear(self):
        return True

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        # calculate dc parameters
        self.calc_dc(x)

        B = self.n1
        C = self.n2
        E = self.n3
        Vbe = self.oppoint['Vbe']
        Vbc = self.oppoint['Vbc']
        gpi = self.oppoint['gpi']
        gmu = self.oppoint['gmu']
        gmf = self.oppoint['gmf']
        gmr = self.oppoint['gmr']
        Ibe = self.oppoint['Ibe']
        Ibc = self.oppoint['Ibc']
        It = self.oppoint['It']

        Ibeeq = Ibe - gpi * Vbe
        Ibceq = Ibc - gmu * Vbc
        Iceeq = It  - gmf * Vbe + gmr * Vbc

        # fill MNA matrices
        A[B][B] = A[B][B] + gmu + gpi
        A[B][C] = A[B][C] - gmu
        A[B][E] = A[B][E] - gpi
        A[C][B] = A[C][B] - gmu + gmf - gmr
        A[C][C] = A[C][C] + gmu + gmr
        A[C][E] = A[C][E] - gmf
        A[E][B] = A[E][B] - gpi - gmf + gmr
        A[E][C] = A[E][C] - gmr
        A[E][E] = A[E][E] + gpi + gmf
        z[B] = z[B] - Ibeeq - Ibceq
        z[C] = z[C] + Ibceq - Iceeq
        z[E] = z[E] + Ibeeq + Iceeq

    def add_ac_stamps(self, A, z, x, iidx, freq):
        B = self.n1
        C = self.n2
        E = self.n3
        S = self.n4
        gpi = self.oppoint['gpi']
        gmu = self.oppoint['gmu']
        gmf = self.oppoint['gmf']
        gmr = self.oppoint['gmr']
        Cbci = self.oppoint['Cbci']
        Cbe = self.oppoint['Cbe']
        Ccs = self.oppoint['Ccs']
        Cbebc = self.oppoint['Cbebc']
        Tf = self.options['Tf']
        Ptf = self.options['Ptf']

        omega = 2. * np.pi * freq
        gmf = gmf * np.exp(-1j * (np.pi / 180. * Ptf) * omega * Tf)
        Ybc = gmu + 1j * omega * Cbci
        Ybe = gpi + 1j * omega * Cbe
        Ycs = 1j * omega * Ccs
        Ybebc = 1j * omega * Cbebc

        # fill admittance matrix
        A[B][B] = A[B][B] + Ybc + Ybe + Ybebc
        A[B][C] = A[B][C] - Ybc - Ybebc
        A[B][E] = A[B][E] - Ybe
        A[B][S] = A[B][S] + 0.
        A[C][B] = A[C][B] + gmf - Ybc - gmr
        A[C][C] = A[C][C] + Ycs + Ybc + gmr
        A[C][E] = A[C][E] - gmf
        A[C][S] = A[C][S] - Ycs
        A[E][B] = A[E][B] + gmr - gmf - Ybe - Ybebc
        A[E][C] = A[E][C] - gmr + Ybebc
        A[E][E] = A[E][E] + gmf + Ybe
        A[E][S] = A[E][S] + 0.
        A[S][B] = A[S][B] + 0.
        A[S][C] = A[S][C] - Ycs
        A[S][E] = A[S][E] + 0.
        A[S][S] = A[S][S] + Ycs

    def calc_oppoint(self, x):
        self.init()
        self.calc_dc(x)

        Cje = self.options['Cje']
        Vje = self.options['Vje']
        Mje = self.options['Mje']
        Cjc = self.options['Cjc']
        Vjc = self.options['Vjc']
        Mjc = self.options['Mjc']
        Cjs = self.options['Cjs']
        Vjs = self.options['Vjs']
        Mjs = self.options['Mjs']
        Fc = self.options['Fc']
        Xcjc = self.options['Xcjc']
        Tf = self.options['Tf']
        Xtf = self.options['Xtf']
        Vtf = self.options['Vtf']
        Itf = self.options['Itf']
        Tr = self.options['Tr']

        Vbe = self.oppoint['Vbe']
        Vbc = self.oppoint['Vbc']
        Vsc = self.oppoint['Vsc']
        gpi = self.oppoint['gpi']
        gmu = self.oppoint['gmu']
        gif = self.oppoint['gif']
        gir = self.oppoint['gir']
        gmf = self.oppoint['gmf']
        gmr = self.oppoint['gmr']
        Ibe = self.oppoint['Ibe']
        Ibc = self.oppoint['Ibc']
        If = self.oppoint['If']
        Qb = self.oppoint['Qb']
        dQb_dVbe = self.oppoint['dQb_dVbe']
        dQb_dVbc = self.oppoint['dQb_dVbc']

        if Vbe <= Fc * Vje:
            Cbedep = Cje * np.power((1. - (Vbe / Vje)), -Mje)
        else:
            Cbedep = Cje / np.power((1. - Fc), Mje) * (1. + Mje * (Vbe - Fc * Vje) / (Vje * (1. - Fc)))

        if Vbc <= Fc * Vjc:
            Cbcdep = Cjc * np.power((1. - (Vbc / Vjc)), -Mjc)
        else:
            Cbcdep = Cjc / np.power((1. - Fc), Mjc) * (1. + Mjc * (Vbc - Fc * Vjc) / (Vjc * (1. - Fc)))

        if Vsc <= 0:
            Cscdep = Cjs * np.power((1. - (Vsc / Vjs)), -Mjs)
        else:
            Cscdep = Cjs * (1. + Mjs * Vsc / Vjs)

        Tff = Tf * (1. + Xtf * np.square(If / (If + Itf)) * np.exp(Vbc / (1.44 * Vtf)))
        dTff_dVbe = Tf * Xtf * 2 * gif * If * Itf * np.exp(Vbc / (1.44 * Vtf))
        dTff_dVbc = (Tf * Xtf / (1.44 * Vtf)) * np.square(If / (If + Itf)) * np.exp(Vbc / (1.44 * Vtf))

        Cbcidep = Xcjc * Cbcdep
        Cbcxdep = (1. - Xcjc) * Cbcdep 
        Cbcdiff = Tr * gir
        Cbediff = (1. / Qb) * (If * dTff_dVbe + Tff * (gif - (If / Qb) * dQb_dVbe))
        Cbebc = (If / Qb) * (dTff_dVbc - (Tff / Qb) * dQb_dVbc)

        self.oppoint['Cbci'] = Cbcidep + Cbcdiff
        self.oppoint['Cbe'] = Cbedep + Cbediff
        self.oppoint['Ccs'] = Cscdep
        self.oppoint['Cbebc'] = Cbebc

    def calc_dc(self, x):
        B   = self.n1
        C   = self.n2
        E   = self.n3
        S   = self.n4
        Vt  = k * self.options['Temp'] / e
        Is  = self.options['Is'] 
        Nf  = self.options['Nf'] 
        Nr  = self.options['Nr'] 
        Ikf = self.options['Ikf']
        Ikr = self.options['Ikr']
        Vaf = self.options['Vaf']
        Var = self.options['Var']
        Ise = self.options['Ise']
        Ne  = self.options['Ne'] 
        Isc = self.options['Isc']
        Nc  = self.options['Nc'] 
        Bf  = self.options['Bf'] 
        Br  = self.options['Br'] 
        Rbm = self.options['Rbm']
        Irb = self.options['Irb']
        Rc  = self.options['Rc'] 
        Re  = self.options['Re'] 
        Rb  = self.options['Rb'] 

        Vb = x[B-1] if B > 0 else 0.0
        Vc = x[C-1] if C > 0 else 0.0
        Ve = x[E-1] if E > 0 else 0.0
        Vs = x[S-1] if S > 0 else 0.0
        Vbe = Vb - Ve
        Vbc = Vb - Vc
        Vsc = Vs - Vc
        
        Vbe, Vbc = self.limit_bjt_voltages(Vbe, Vbc, Vt)

        If = Is * np.expm1(Vbe / (Nf * Vt))

        Ibei = If / Bf
        Iben = Ise * np.expm1(Vbe / (Ne * Vt))
        Ibe  = Ibei + Iben

        gbei = Is / (Nf * Vt * Bf) * np.exp(Vbe / (Nf * Vt))
        gben = Ise / (Ne * Vt) * np.exp(Vbe / (Ne * Vt))
        gpi  = gbei + gben

        Ir = Is * np.expm1(Vbc / (Nr * Vt))

        Ibci = Ir / Br
        Ibcn = Isc * np.expm1(Vbc / (Nc * Vt))
        Ibc  = Ibci + Ibcn

        gbci = Is / (Nr * Vt * Br) * np.exp(Vbc / (Nr * Vt))
        gbcn = Isc / (Nc * Vt) * np.exp(Vbc / (Nc * Vt))
        gmu  = gbci + gbcn

        Q2 = (If / Ikf) + (Ir / Ikr)
        Q1 = 1. / (1. - (Vbc / Vaf) - (Vbe / Var))
        Qb = (Q1 / 2.) * (1. + np.sqrt(1. + 4. * Q2))

        It = (If - Ir) / Qb

        gif = gbei * Bf
        gir = gbci * Br

        dQb_dVbe = Q1 * ((Qb / Var) + (gif / (Ikf * np.sqrt(1. + 4. * Q2))))
        dQb_dVbc = Q1 * ((Qb / Vaf) + (gir / (Ikr * np.sqrt(1. + 4. * Q2))))

        gmf = (1. / Qb) * (+ gif - It * dQb_dVbe)
        gmr = (1. / Qb) * (- gir - It * dQb_dVbc)

        # save dc point parameters
        self.oppoint['Vb'] = Vb
        self.oppoint['Vc'] = Vc
        self.oppoint['Ve'] = Ve
        self.oppoint['Vs'] = Vs
        self.oppoint['Vbe'] = Vbe
        self.oppoint['Vbc'] = Vbc
        self.oppoint['Vsc'] = Vsc
        self.oppoint['If'] = If
        self.oppoint['Ir'] = Ir
        self.oppoint['Ib'] = Ibe + Ibc
        self.oppoint['Ic'] = It - Ibc
        self.oppoint['Ie'] = It + Ibe
        self.oppoint['Ibe'] = Ibe
        self.oppoint['Ibc'] = Ibc
        self.oppoint['It'] = It
        self.oppoint['gpi'] = gpi
        self.oppoint['gmu'] = gmu
        self.oppoint['gif'] = gif
        self.oppoint['gir'] = gir
        self.oppoint['gmf'] = gmf
        self.oppoint['gmr'] = gmr
        self.oppoint['Qb'] = Qb
        self.oppoint['dQb_dVbe'] = dQb_dVbe
        self.oppoint['dQb_dVbc'] = dQb_dVbc

    def limit_bjt_voltages(self, Vbe, Vbc, Vt):
        Is = self.options['Is']
        Nf = self.options['Nf']
        Nr = self.options['Nr']

        # limiting algorithm to avoid overflow
        Vbecrit = Nf * Vt * np.log(Nf * Vt / (np.sqrt(2.) * Is))
        Vbccrit = Nr * Vt * np.log(Nr * Vt / (np.sqrt(2.) * Is))
        if (Vbe > Vbecrit and Vbe > 0.):
            Vbe = self.Vbeold + Nf * Vt * np.log1p((Vbe - self.Vbeold) / (Nf * Vt))
        if (Vbc > Vbccrit and Vbc > 0.):
            Vbc = self.Vbcold + Nr * Vt * np.log1p((Vbc - self.Vbcold) / (Nr * Vt))
        self.Vbeold = Vbe
        self.Vbcold = Vbc

        return Vbe, Vbc

    def __str__(self):
        return 'BJT: {}\nNodes BCE nodes = {}, {}, {}\n'.format(self.name, self.n1, self.n2, self.n3)
