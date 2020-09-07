import numpy as np

class TransientVoltageSource():
    
    def __init__(self, name, n1, n2, dc=0, tstart=0, tstop=1e3, vtype=None):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.dc = float(dc)
        self.tstart = float(tstart)
        self.tstop = float(tstop)
        self.vtype = vtype

    def get_num_vsources(self, analysis):
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
        z[iidx] = 0 # 0V

    def add_ac_stamps(self, A, z, x, iidx, freq):
        A[self.n1][iidx] = +1.0 
        A[self.n2][iidx] = -1.0 
        A[iidx][self.n1] = +1.0 
        A[iidx][self.n2] = -1.0
        z[iidx] = 0 # 0V

    def add_tran_stamps(self, A, z, x, iidx, t, tstep):
        if self.vtype == 'step':
            A[self.n1][iidx] = +1.0 
            A[self.n2][iidx] = -1.0 
            A[iidx][self.n1] = +1.0 
            A[iidx][self.n2] = -1.0
            if t < self.tstart or t > self.tstop:
                z[iidx] = 0
            else:
                z[iidx] = self.dc

    def __str__(self):
        return 'Voltage Source: {}\nNodes: {} -> {}\nVdc = {}\nVac = {}\nPhase = {}'.format(self.name,
                                                                                            self.n1,
                                                                                            self.n2,
                                                                                            self.dc,
                                                                                            self.ac,
                                                                                            np.degrees(self.phase))

