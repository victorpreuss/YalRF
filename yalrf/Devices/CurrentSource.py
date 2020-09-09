import numpy as np

class CurrentSource():
    
    def __init__(self, name, n1, n2, dc=0, ac=0, phase=0, itype=None):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.itype = itype
        self.dc = float(dc)
        self.ac = float(ac)
        self.phase = np.radians(float(phase))

    def get_num_vsources(self, analysis):
        return 0

    def is_nonlinear(self):
        return False

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        if self.itype != 'ac': # ac current source is open (0A) at dc
            z[self.n1] = z[self.n1] + self.dc
            z[self.n2] = z[self.n2] - self.dc

    def add_ac_stamps(self, A, z, x, iidx, freq):
        if self.itype != 'dc': # dc current source is open (0A) at ac
            z[self.n1] = z[self.n1] + self.ac * np.exp(1j * self.phase)
            z[self.n2] = z[self.n2] - self.ac * np.exp(1j * self.phase)

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        self.add_dc_stamps(A, z, x, iidx)

    def __str__(self):
        return 'Current Source: {}\nNodes = {} -> {}\nIdc = {}\nIac = {}\nPhase = {}'.format(self.name,
                                                                                             self.n1,
                                                                                             self.n2,
                                                                                             self.dc,
                                                                                             self.ac,
                                                                                             np.degrees(self.phase))

