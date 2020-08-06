class VoltageSource():
    
    def __init__(self, name, n1, n2, vtype=None, dc=None, mag=None, phase=None):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.vtype = vtype
        self.dc = dc
        self.mag = mag
        self.phase = phase

    def get_vsource(self):
        return 1

    def is_linear(self):
        return True

    def add_dc_stamps(self, A, z, x, iidx):
        A[self.n1][iidx] = +1.0 
        A[self.n2][iidx] = -1.0 
        A[iidx][self.n1] = +1.0 
        A[iidx][self.n2] = -1.0
        z[iidx] = self.V()

    def V(self):
        if self.vtype.lower() == 'dc':
            return self.dc
        else:
            print('WARNING: unknown voltage source')
            return 0

    def __str__(self):
        out = ''
        if self.vtype.lower() == 'dc':
            out = 'Voltage Source: {}\nNodes = {} -> {}\nVdc = {}\n'.format(self.name, self.n1, self.n2, self.dc)
        return out

