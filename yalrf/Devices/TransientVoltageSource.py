import numpy as np

class TransientVoltageSource():
    
    def __init__(self, name, n1, n2, vtype='sine',
                                     dc=0, ac=0, freq=0, phase=0,
                                     v1=0, v2=0, tstart=0, tstop=1e3, trise=1e-12, tfall=1e-12):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.vtype = vtype

        # sine
        self.dc = float(dc)
        self.ac = float(ac)
        self.freq = float(freq)
        self.phase = np.radians(phase)

        # pulse
        self.v1 = float(v1)
        self.v2 = float(v2)
        self.tstart = float(tstart)
        self.tstop = float(tstop)
        self.tstart = float(tstart)
        self.trise = float(trise)
        self.tfall = float(tfall)

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

        if self.vtype == 'sine':
            z[iidx] = self.dc + self.ac * np.sin(self.phase)
        elif self.vtype == 'pulse':
            z[iidx] = self.v1
        else:
            z[iidx] = 0.

    def add_ac_stamps(self, A, z, x, iidx, freq):
        A[self.n1][iidx] = +1.0
        A[self.n2][iidx] = -1.0
        A[iidx][self.n1] = +1.0
        A[iidx][self.n2] = -1.0

        if self.vtype == 'sine':
            z[iidx] = self.ac
        else:
            z[iidx] = 0.

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        A[self.n1][iidx] = +1.0
        A[self.n2][iidx] = -1.0
        A[iidx][self.n1] = +1.0
        A[iidx][self.n2] = -1.0

        if self.vtype == 'sine':
            z[iidx] = self.dc + self.ac * np.sin(2. * np.pi * self.freq * t + self.phase)
        elif self.vtype == 'pulse':
            if t < self.tstart:
                z[iidx] = self.v1
            elif t < self.tstart + self.trise:
                z[iidx] = self.v1 + ((self.v2 - self.v1) / self.trise) * (t - self.tstart)
            elif t < self.tstop:
                z[iidx] = self.v2
            elif t < self.tstop + self.tfall:
                z[iidx] = self.v2 + ((self.v1 - self.v2) / self.tfall) * (t - self.tstop)
            else:
                z[iidx] = self.v1

    def __str__(self):
        if self.vtype == 'sine':
            return 'Sinusoidal Source: {}\nNodes: {} -> {}\nVdc = {}\nVac = {}\nFreq = {}\nPhase = {}'.format(self.name,
                                                                                                              self.n1,
                                                                                                              self.n2,
                                                                                                              self.dc,
                                                                                                              self.ac,
                                                                                                              self.freq,
                                                                                                              np.degrees(self.phase))
        elif self.vtype == 'pulse':
            return 'Pulse Source: {}\nNodes: {} -> {}\nV1 = {}\nV2 = {}\nTstart = {}\nTstop = {}'.format(self.name,
                                                                                                         self.n1,
                                                                                                         self.n2,
                                                                                                         self.v1,
                                                                                                         self.v2,
                                                                                                         self.tstart,
                                                                                                         self.tstop)
        else:
            return 'None'

