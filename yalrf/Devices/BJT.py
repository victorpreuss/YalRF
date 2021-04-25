from scipy.constants import k, e
import numpy as np

# BJT Spice Gummel-Poon model options
options = {}
options['Temp'] = 300.0 # device temperature
options['Is']  = 1e-15  # saturation current
options['Nf']  = 1.0    # forward emission coefficient
options['Nr']  = 1.0    # reverse emission coefficient
options['Ikf'] = 1e9    # high current corner for forward beta
options['Ikr'] = 1e9    # high current corner for reverse beta
options['Vaf'] = 1e12   # forward Early Voltage
options['Var'] = 1e12   # reverse Early Voltage
options['Ise'] = 0.0    # base-emitter leakage saturation current
options['Ne']  = 1.5    # base-emitter leakage emission coefficient
options['Isc'] = 0.0    # base-collector leakage saturation current
options['Nc']  = 2.0    # base-collector leakage emission coefficient
options['Bf']  = 100.0  # forward beta
options['Br']  = 1.0    # reverse beta
options['Rbm'] = 0.0    # minimum base resistance for high-currents
options['Irb'] = 1e9    # current for base resistance midpoint
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
options['Vtf'] = 1e12   # coefficient of base-collector voltage dependence for Tf
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
#       excess phase (for transient?)
#       temperature dependence
#       Cbcxdep is currently unused
#       noise
class BJT():
    
    def __init__(self, name, n1, n2, n3, n4=0):
        self.name = name
        self.n1   = n1 # base (B)
        self.n2   = n2 # collector (C)
        self.n3   = n3 # emitter (E)
        self.n4   = n4 # substrate (S)

        self.options = options.copy() # bjt options
        self.oppoint = {}

        # this options vector holds the corrected values
        # for the BJT according to area and temperature
        self.adjusted_options = {}

        # for the limiting scheme
        self.Vbeold = 0.
        self.Vbcold = 0.

        self.VbeoldHB = None
        self.VbcoldHB = None

        self.Ib = []
        self.Ic = []

    def get_num_vsources(self):
        return 0

    def is_nonlinear(self):
        return True

    def get_idc(self, x):
        return self.Ib[0], self.Ic[0], self.Ib[0] + self.Ic[0]

    def get_itran(self, x):
        Ib = np.array(self.Ib)
        Ic = np.array(self.Ic)
        Ie = Ib + Ic
        return Ib, Ic, Ie

    def init(self):
        # clear model internal variables
        self.oppoint = {}

        self.Vbeold = 0.
        self.Vbcold = 0.
        self.Ib = []
        self.Ic = []
        
        # area and temperature dependent adjustments
        A = self.options['Area']

        self.adjusted_options['Is'] = self.options['Is'] * A
        self.adjusted_options['Ise'] = self.options['Ise'] * A
        self.adjusted_options['Isc'] = self.options['Isc'] * A
        self.adjusted_options['Ikf'] = self.options['Ikf'] * A
        self.adjusted_options['Ikr'] = self.options['Ikr'] * A
        self.adjusted_options['Irb'] = self.options['Irb'] * A
        self.adjusted_options['Itf'] = self.options['Itf'] * A

        self.adjusted_options['Cje'] = self.options['Cje'] * A
        self.adjusted_options['Cjs'] = self.options['Cjs'] * A
        self.adjusted_options['Cjc'] = self.options['Cjc'] * A

        self.adjusted_options['Rb'] = self.options['Rb'] / A
        self.adjusted_options['Rbm'] = self.options['Rbm'] / A
        self.adjusted_options['Rc'] = self.options['Rc'] / A
        self.adjusted_options['Re'] = self.options['Re'] / A

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
        Iceeq = It  - gmf * Vbe - gmr * Vbc

        # fill MNA matrices
        A[B][B] = A[B][B] + gmu + gpi
        A[B][C] = A[B][C] - gmu
        A[B][E] = A[B][E] - gpi
        A[C][B] = A[C][B] - gmu + gmf + gmr
        A[C][C] = A[C][C] + gmu - gmr
        A[C][E] = A[C][E] - gmf
        A[E][B] = A[E][B] - gpi - gmf - gmr
        A[E][C] = A[E][C] + gmr
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

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        self.add_dc_stamps(A, z, x, iidx)

    def add_hb_stamps(self, v, i, g, k):
        B = self.n1 
        C = self.n2 
        E = self.n3 

        # calculate currents at the BJT
        Vb = v[B-1,k] if B > 0 else 0
        Vc = v[C-1,k] if C > 0 else 0
        Ve = v[E-1,k] if E > 0 else 0
        Vs = 0

        Ib, Ic, Ie = self.get_i(Vb, Vc, Ve, Vs)
        gmu, gpi, gmf, gmr = self.get_g(Vb, Vc, Ve, Vs)

        g[B,B,k] = g[B,B,k] + gmu + gpi
        g[B,C,k] = g[B,C,k] - gmu
        g[B,E,k] = g[B,E,k] - gpi
        g[C,B,k] = g[C,B,k] - gmu + gmf + gmr
        g[C,C,k] = g[C,C,k] + gmu - gmr
        g[C,E,k] = g[C,E,k] - gmf
        g[E,B,k] = g[E,B,k] - gpi - gmf - gmr
        g[E,C,k] = g[E,C,k] + gmr
        g[E,E,k] = g[E,E,k] + gpi + gmf

        i[B,k] = i[B,k] + Ib
        i[C,k] = i[C,k] + Ic
        i[E,k] = i[E,k] - Ie

    def save_oppoint(self):
        Ib = self.oppoint['Ib']
        Ic = self.oppoint['Ic']
        self.Ib.append(Ib)
        self.Ic.append(Ic)

    def save_tran(self, x, tstep):
        Ib = self.oppoint['Ib']
        Ic = self.oppoint['Ic']
        self.Ib.append(Ib)
        self.Ic.append(Ic)

    def calc_oppoint(self, x):
        self.calc_dc(x)

        Cje = self.adjusted_options['Cje']
        Vje = self.options['Vje']
        Mje = self.options['Mje']
        Cjc = self.adjusted_options['Cjc']
        Vjc = self.options['Vjc']
        Mjc = self.options['Mjc']
        Cjs = self.adjusted_options['Cjs']
        Vjs = self.options['Vjs']
        Mjs = self.options['Mjs']
        Fc = self.options['Fc']
        Xcjc = self.options['Xcjc']
        Tf = self.options['Tf']
        Xtf = self.options['Xtf']
        Vtf = self.options['Vtf']
        Itf = self.adjusted_options['Itf']
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
        Is  = self.adjusted_options['Is'] 
        Nf  = self.options['Nf'] 
        Nr  = self.options['Nr'] 
        Ikf = self.adjusted_options['Ikf']
        Ikr = self.adjusted_options['Ikr']
        Vaf = self.options['Vaf']
        Var = self.options['Var']
        Ise = self.adjusted_options['Ise']
        Ne  = self.options['Ne'] 
        Isc = self.adjusted_options['Isc']
        Nc  = self.options['Nc'] 
        Bf  = self.options['Bf'] 
        Br  = self.options['Br'] 
        Rbm = self.adjusted_options['Rbm']
        Irb = self.adjusted_options['Irb']
        Rc  = self.adjusted_options['Rc'] 
        Re  = self.adjusted_options['Re'] 
        Rb  = self.adjusted_options['Rb'] 

        Vb = x[B-1,0] if B > 0 else 0.
        Vc = x[C-1,0] if C > 0 else 0.
        Ve = x[E-1,0] if E > 0 else 0.
        Vs = x[S-1,0] if S > 0 else 0.
        Vbe = Vb - Ve
        Vbc = Vb - Vc
        Vsc = Vs - Vc

        gmin = 1e-12
        
        Vbe, Vbc = self.limit_bjt_voltages(Vbe, Vbc, Vt)

        If = Is * (np.exp(Vbe / (Nf * Vt)) - 1.)

        Ibei = If / Bf
        Iben = Ise * (np.exp(Vbe / (Ne * Vt)) - 1.) + gmin * Vbe
        Ibe  = Ibei + Iben

        gbei = Is / (Nf * Vt * Bf) * np.exp(Vbe / (Nf * Vt))
        gben = Ise / (Ne * Vt) * np.exp(Vbe / (Ne * Vt)) + gmin
        gpi  = gbei + gben

        Ir = Is * (np.exp(Vbc / (Nr * Vt)) - 1.)

        Ibci = Ir / Br
        Ibcn = Isc * (np.exp(Vbc / (Nc * Vt)) - 1.) + gmin * Vbc
        Ibc  = Ibci + Ibcn

        gbci = Is / (Nr * Vt * Br) * np.exp(Vbc / (Nr * Vt))
        gbcn = Isc / (Nc * Vt) * np.exp(Vbc / (Nc * Vt)) + gmin
        gmu  = gbci + gbcn

        Q1 = 1. / (1. - (Vbc / Vaf) - (Vbe / Var))
        Q2 = (If / Ikf) + (Ir / Ikr)
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

    def check_vlimit(self, x, vabstol):
        B   = self.n1
        C   = self.n2
        E   = self.n3
        
        Vt  = k * self.options['Temp'] / e
        Vb = x[B-1,0] if B > 0 else 0.
        Vc = x[C-1,0] if C > 0 else 0.
        Ve = x[E-1,0] if E > 0 else 0.

        Vbe = Vb - Ve
        Vbc = Vb - Vc

        Vbelim, Vbclim = self.limit_bjt_voltages(Vbe, Vbc, Vt)

        if (Vbe - Vbelim > vabstol) or (Vbc - Vbclim > vabstol):
            return False

        return True

    def limit_bjt_voltages(self, Vbe, Vbc, Vt):
        Is = self.options['Is']
        Nf = self.options['Nf']
        Nr = self.options['Nr']

        # limiting algorithm to avoid overflow
        # Vbecrit = Nf * Vt * np.log(Nf * Vt / (np.sqrt(2.) * Is))
        # Vbccrit = Nr * Vt * np.log(Nr * Vt / (np.sqrt(2.) * Is))
        # if (Vbe > Vbecrit and Vbe > 0.):
        #     Vbe = self.Vbeold + Nf * Vt * np.log1p((Vbe - self.Vbeold) / (Nf * Vt))
        # if (Vbc > Vbccrit and Vbc > 0. and ((Vbc - self.Vbcold) / (Nr * Vt)) > -1.):
        #     Vbc = self.Vbcold + Nr * Vt * np.log1p((Vbc - self.Vbcold) / (Nr * Vt))

        Vbe = self.Vbeold + 10. * Nf * Vt * np.tanh((Vbe - self.Vbeold) / (10. * Nf * Vt))
        Vbc = self.Vbcold + 10. * Nr * Vt * np.tanh((Vbc - self.Vbcold) / (10. * Nr * Vt))

        self.Vbeold = Vbe
        self.Vbcold = Vbc

        return Vbe, Vbc

    # TODO: get_hb_params() is a temporary implementations to test
    #       the harmonic balance algorithm. Eventually the complete
    #       BJT model should be used.
    def get_hb_params(self, Vb, Vc, Ve, Vs, s):
        Vt  = k * self.options['Temp'] / e
        Is  = self.adjusted_options['Is'] 
        Nf  = self.options['Nf'] 
        Nr  = self.options['Nr'] 
        Ikf = self.adjusted_options['Ikf']
        Ikr = self.adjusted_options['Ikr']
        Vaf = self.options['Vaf']
        Var = self.options['Var']
        Ise = self.adjusted_options['Ise']
        Ne  = self.options['Ne'] 
        Isc = self.adjusted_options['Isc']
        Nc  = self.options['Nc'] 
        Bf  = self.options['Bf'] 
        Br  = self.options['Br'] 
        Cje = self.adjusted_options['Cje']
        Vje = self.options['Vje']
        Mje = self.options['Mje']
        Cjc = self.adjusted_options['Cjc']
        Vjc = self.options['Vjc']
        Mjc = self.options['Mjc']
        Cjs = self.adjusted_options['Cjs']
        Vjs = self.options['Vjs']
        Mjs = self.options['Mjs']
        Fc = self.options['Fc']
        Xcjc = self.options['Xcjc']
        Tf = self.options['Tf']
        Xtf = self.options['Xtf']
        Vtf = self.options['Vtf']
        Itf = self.adjusted_options['Itf']
        Tr = self.options['Tr']

        # TODO: improve this later
        Vbe = Vb - Ve + 1e-6
        Vbc = Vb - Vc + 1e-6
        Vsc = Vs - Vc + 1e-6

        # Vbe = self.VbeoldHB[s] + 10. * Nf * Vt * np.tanh((Vbe - self.VbeoldHB[s]) / (10. * Nf * Vt))
        # Vbc = self.VbcoldHB[s] + 10. * Nr * Vt * np.tanh((Vbc - self.VbcoldHB[s]) / (10. * Nr * Vt))

        self.VbeoldHB[s] = Vbe
        self.VbcoldHB[s] = Vbc

        # print(Vbe)
        # print(Vbc)

        gmin = 1e-12
        cmin = 1e-18
        
        If = Is * (exp_lim(Vbe / (Nf * Vt)) - 1.)

        Ibei = If / Bf
        Iben = Ise * (exp_lim(Vbe / (Ne * Vt)) - 1.)
        Ibe  = Ibei + Iben + gmin * Vbe

        Ir = Is * (exp_lim(Vbc / (Nr * Vt)) - 1.)

        Ibci = Ir / Br
        Ibcn = Isc * (exp_lim(Vbc / (Nc * Vt)) - 1.)
        Ibc  = Ibci + Ibcn + gmin * Vbc

        Q1 = 1. / (1. - (Vbc / Vaf) - (Vbe / Var))
        Q2 = (If / Ikf) + (Ir / Ikr)
        Qb = (Q1 / 2.) * (1. + np.sqrt(1. + 4. * Q2))

        It = (If - Ir) / Qb

        gbei = Is / (Nf * Vt * Bf) * exp_lim(Vbe / (Nf * Vt))
        gben = Ise / (Ne * Vt) * exp_lim(Vbe / (Ne * Vt))
        gpi  = gbei + gben + gmin

        gbci = Is / (Nr * Vt * Br) * exp_lim(Vbc / (Nr * Vt))
        gbcn = Isc / (Nc * Vt) * exp_lim(Vbc / (Nc * Vt))
        gmu  = gbci + gbcn + gmin

        gif = gbei * Bf
        gir = gbci * Br

        dQb_dVbe = Q1 * ((Qb / Var) + (gif / (Ikf * np.sqrt(1. + 4. * Q2))))
        dQb_dVbc = Q1 * ((Qb / Vaf) + (gir / (Ikr * np.sqrt(1. + 4. * Q2))))

        gmf = (1. / Qb) * (+ gif - It * dQb_dVbe)
        gmr = (1. / Qb) * (- gir - It * dQb_dVbc)

        Cbedep = pn_capacitance(Vbe, Cje, Vje, Mje, Fc)
        Qbe = pn_charge(Vbe, Cje, Vje, Mje, Fc)

        Cbcdep = pn_capacitance(Vbc, Cjc, Vjc, Mjc, Fc)
        Qbc = pn_charge(Vbc, Cjc, Vjc, Mjc, Fc)

        if Vsc <= 0:
            Cscdep = Cjs * np.power((1. - (Vsc / Vjs)), -Mjs)
            Qsc = Cjs * Vjs / (1 - Mjs) * (1 - np.power(1. - Vsc / Vjs, 1. - Mjs))
        else:
            Cscdep = Cjs * (1. + Mjs * Vsc / Vjs)
            Qsc = Cjs * Vsc / (1 + (Mjs * Vsc) / (2. * Vjs))

        Tff = Tf * (1. + Xtf * np.square(If / (If + Itf)) * exp_lim(Vbc / (1.44 * Vtf)))
        dTff_dVbe = Tf * Xtf * 2 * gif * If * Itf * exp_lim(Vbc / (1.44 * Vtf))
        dTff_dVbc = (Tf * Xtf / (1.44 * Vtf)) * np.square(If / (If + Itf)) * exp_lim(Vbc / (1.44 * Vtf))

        Cbcidep = Xcjc * Cbcdep
        Cbcxdep = (1. - Xcjc) * Cbcdep 
        Cbcdiff = Tr * gir
        Cbediff = (1. / Qb) * (If * dTff_dVbe + Tff * (gif - (If / Qb) * dQb_dVbe))
        Cbebc = (If / Qb) * (dTff_dVbc - (Tff / Qb) * dQb_dVbc) + cmin

        Ib = Ibe + Ibc
        Ic = It - Ibc
        Ie = Ib + Ic

        if Xcjc != 1.:
            print('WARNING: external base-collector capacitance is currently unsupported')

        Cbc = Cbcidep + Cbcdiff + cmin
        Qbc += Ir * Tr + cmin * Vbc

        Cbe = Cbedep + Cbediff + cmin
        Qbe += If * Tff / Qb + cmin * Vbe

        Csc = Cscdep + cmin
        Qsc += cmin * Vsc

        return Ib, Ic, Ie, Qbe, Qbc, Qsc, gmu, gpi, gmf, gmr, Cbc, Cbe, Cbebc, Csc

    def __str__(self):
        return 'BJT: {}\nNodes BCE nodes = {}, {}, {}\n'.format(self.name, self.n1, self.n2, self.n3)

# limit the maximum derivative of the exponential function
def exp_lim(x):
    return np.exp(x) if x < 200. else np.exp(200.) + np.exp(200.) * (x - 200.)

def pn_capacitance(Vpn, Cj, Vj, Mj, Fc):
    if Vpn <= Fc * Vj:
        C = Cj * np.power((1. - (Vpn / Vj)), -Mj)
    else:
        C = Cj / np.power((1. - Fc), Mj) * (1. + Mj * (Vpn - Fc * Vj) / (Vj * (1. - Fc)))

    return C

def pn_charge(Vpn, Cj, Vj, Mj, Fc):
    if Vpn <= Fc * Vj:
        Q = Cj * Vj / (1. - Mj) * (1. - np.power((1. - Vpn / Vj), (1. - Mj)))
    else:
        X = (1. - np.power((1. - Fc), (1. - Mj))) / (1. - Mj) + \
            (1. - Fc * (1. + Mj)) / np.power((1. - Fc), (1. + Mj)) * (Vpn / Vj - Fc) + \
            Mj / (2. * np.power((1. - Fc), (1. + Mj))) * (np.square(Vpn / Vj) - np.square(Fc))
        Q = Cj * Vj * X

    return Q

"""

    # LEGACY CODE TO USE SECANT METHOD FOR dIdV

    # USING SECANT METHOD 
    Ibn, Icn, Ien = dev.get_i(Vb * 1.01, Vc, Ve, Vs)

    gbb = (Ibn - Ib) / (0.01 * Vb)
    gcb = (Icn - Ic) / (0.01 * Vb)
    geb = (Ien - Ie) / (0.01 * Vb)

    gt[B,B,s] = gt[B,B,s] + gbb
    gt[C,B,s] = gt[C,B,s] + gcb
    gt[E,B,s] = gt[E,B,s] + geb

    Ibn, Icn, Ien = dev.get_i(Vb, Vc * 1.01, Ve, Vs)

    gcc = (Icn - Ic) / (0.01 * Vc)
    gbc = (Ibn - Ib) / (0.01 * Vc)
    gec = (Ien - Ie) / (0.01 * Vc)

    gt[C,C,s] = gt[C,C,s] + gcc
    gt[B,C,s] = gt[B,C,s] + gbc
    gt[E,C,s] = gt[E,C,s] + gec

    Ibn, Icn, Ien = dev.get_i(Vb, Vc, Ve * 1.01, Vs)

    gee = (Ien - Ie) / (0.01 * Ve)
    gbe = (Ibn - Ib) / (0.01 * Ve)
    gce = (Icn - Ic) / (0.01 * Ve)

    gt[E,E,s] = gt[E,E,s] + gee
    gt[B,E,s] = gt[B,E,s] + gbe
    gt[C,E,s] = gt[C,E,s] + gce

    it[B,s] = it[B,s] + Ib
    it[C,s] = it[C,s] + Ic
    it[E,s] = it[E,s] - Ie

"""
