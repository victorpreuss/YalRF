import numpy as np

class VoltageSource():
    
    def __init__(self, name, n1, n2, dc=0, ac=0, phase=0, vtype=None):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.vtype = vtype
        self.dc = float(dc)
        self.ac = float(ac)
        self.phase = np.radians(float(phase))

    def get_num_vsources(self):
        return 1

    def is_nonlinear(self):
        return False

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        A[self.n1][iidx] = +1.0
        A[self.n2][iidx] = -1.0
        A[iidx][self.n1] = +1.0
        A[iidx][self.n2] = -1.0
        z[iidx] = self.dc

        # ac voltage source is a short (0V) at dc
        if self.vtype == 'ac':
            z[iidx] = 0        

    def add_ac_stamps(self, A, z, x, iidx, freq):
        A[self.n1][iidx] = +1.0 
        A[self.n2][iidx] = -1.0 
        A[iidx][self.n1] = +1.0 
        A[iidx][self.n2] = -1.0
        z[iidx] = self.ac * np.exp(1j * self.phase)

        # dc voltage source is a short (0V) at ac
        if self.vtype == 'dc':
            z[iidx] = 0

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        self.add_dc_stamps(A, z, x, iidx)

    def __str__(self):
        return 'Voltage Source: {}\nNodes: {} -> {}\nVdc = {}\nVac = {}\nPhase = {}'.format(self.name,
                                                                                            self.n1,
                                                                                            self.n2,
                                                                                            self.dc,
                                                                                            self.ac,
                                                                                            np.degrees(self.phase))

