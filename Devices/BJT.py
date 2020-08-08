from scipy.constants import k, e
import numpy as np

# Disclaimer: this implementation was entirely done using
#             the QucsStudio documentation as reference

# BJT Spice Gummel-Poon model options
options = dict()
options['Temp'] = 300.0 # model temperature
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
options['Rbm'] = 0.0
options['Irb'] = 0.0
options['Rc']  = 0.0
options['Re']  = 0.0
options['Rb']  = 0.0

class BJT():
    
    def __init__(self, name, n1, n2, n3):
        self.name = name
        self.n1   = n1 # base (B)
        self.n2   = n2 # collector (C)
        self.n3   = n3 # emitter (E)

        self.options = options.copy() # diode options
        
        # useful for the diodes limiting scheme
        self.Vbeold = 0
        self.Vbcold = 0

    def get_num_vsources(self, analysis):
        return 0

    def is_linear(self):
        return False

    def add_dc_stamps(self, A, z, x, iidx):

        B   = self.n1
        C   = self.n2
        E   = self.n3
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

        Vb = x[B-1] if B else 0.0
        Vc = x[C-1] if C else 0.0
        Ve = x[E-1] if E else 0.0

        Vbe = Vb - Ve
        Vbc = Vb - Vc
        
        Vbecrit = Nf * Vt * np.log(Nf * Vt / (np.sqrt(2) * Is))
        Vbccrit = Nr * Vt * np.log(Nr * Vt / (np.sqrt(2) * Is))
        
        # limiting algorithm to avoid overflow
        if (Vbe > Vbecrit and Vbe > 0):
            Vbe = self.Vbeold + Nf * Vt * np.log(1 + (Vbe - self.Vbeold) / (Nf * Vt))
        if (Vbc > Vbccrit and Vbc > 0):
            Vbc = self.Vbcold + Nr * Vt * np.log(1 + (Vbc - self.Vbcold) / (Nr * Vt))
        self.Vbeold = Vbe
        self.Vbcold = Vbc

        If = Is * (np.exp(Vbe / (Nf * Vt)) - 1)

        Ibei = If / Bf
        Iben = Ise * (np.exp(Vbe / (Ne * Vt)) - 1)
        Ibe  = Ibei + Iben

        gbei = Is / (Nf * Vt * Bf) * np.exp(Vbe / (Nf * Vt))
        gben = Ise / (Ne * Vt) * np.exp(Vbe / (Ne * Vt))
        gpi  = gbei + gben

        Ir = Is * (np.exp(Vbc / (Nr * Vt)) - 1)

        Ibci = Ir / Br
        Ibcn = Isc * (np.exp(Vbc / (Nc * Vt)) - 1)
        Ibc  = Ibci + Ibcn

        gbci = Is / (Nr * Vt * Br) * np.exp(Vbc / (Nr * Vt))
        gbcn = Isc / (Nc * Vt) * np.exp(Vbc / (Nc * Vt))
        gmu  = gbci + gbcn

        Q2 = (If / Ikf) + (Ir / Ikr)
        Q1 = 1 / (1 - (Vbc / Vaf) - (Vbe / Var))
        Qb = (Q1 / 2) * (1 + np.sqrt(1 + 4 * Q2))

        It = (If - Ir) / Qb

        gif = gbei * Bf
        gir = gbci * Br

        dQb_dVbe = Q1 * ((Qb / Var) + (gif / (Ikf * np.sqrt(1 + 4 * Q2))))
        dQb_dVbc = Q1 * ((Qb / Vaf) + (gir / (Ikr * np.sqrt(1 + 4 * Q2))))

        gmf = (1 / Qb) * (+ gif - It * dQb_dVbe)
        gmr = (1 / Qb) * (- gir - It * dQb_dVbc)

        Ibeeq = Ibe - gpi * Vbe
        Ibceq = Ibc - gmu * Vbc
        Iceeq = It  - gmf * Vbe + gmr * Vbc
        
        # print('Vb = {} Vc = {} Ve = {}'.format(Vb, Vc, Ve))
        # print('Vbecrit = {} Vbccrit = {}'.format(Vbecrit, Vbccrit))
        # print('Vbe = {} Vbc = {}'.format(Vbe, Vbc))
        # print('Q2 = {} Q1 = {} Qb = {}'.format(Q2, Q1, Qb))
        # print('It = {} gif = {} gir = {}'.format(It, gif, gir))
        # print('dQb_dVbe = {} dQb_dVbc = {}'.format(dQb_dVbe, dQb_dVbc))
        # print('gpi = {} gmu = {}'.format(gpi, gmu))
        # print('If = {} Ir = {}'.format(If, Ir))
        # print('Ibeeq = {} Ibceq = {} Iceeq = {}'.format(Ibeeq, Ibceq, Iceeq))
        # print()

        # TODO: implement base resistance (bias dependent), the
        #       DC resistances Rb, Rc and Re and the excess phase

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

    def __str__(self):
        return 'BJT: {}\nNodes BCE nodes = {}, {}, {}\n'.format(self.name, self.n1, self.n2, self.n3)
