import numpy as np

class CurrentSource():
    
    def __init__(self, name, n1, n2, itype=None, dc=0, ac=0, phase=0):
        self.name = name
        self.n1    = n1
        self.n2    = n2
        self.itype = itype
        self.dc    = float(dc)
        self.ac    = float(ac)
        self.phase = np.radians(float(phase))

    def get_num_vsources(self, analysis):
        return 0

    def is_linear(self):
        return True

    def add_dc_stamps(self, A, z, x, iidx):
        z[self.n1] = z[self.n1] + self.dc
        z[self.n2] = z[self.n2] - self.dc

    def add_ac_stamps(self, A, z, x, iidx, freq):
        z[self.n1] = z[self.n1] + self.ac * np.exp(1j * self.phase)
        z[self.n2] = z[self.n2] - self.ac * np.exp(1j * self.phase)

    def I(self):
        if self.itype.lower() == 'dc':
            return self.dc
        else:
            print('WARNING: unknown current source')
            return 0

    def __str__(self):
        out = ''
        if self.itype.lower() == 'dc':
            out = 'Current Source: {}\nNodes = {} -> {}\nIdc = {}\n'.format(self.name, self.n1, self.n2, self.dc)
        return out

