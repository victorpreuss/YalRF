class CurrentSource():
    
    def __init__(self, name, n1, n2, itype=None, dc=None, mag=None, phase=None):
        self.name = name
        self.n1    = n1
        self.n2    = n2
        self.itype = itype
        self.dc    = dc
        self.mag   = mag
        self.phase = phase

    def get_vsource(self):
        return 0

    def is_linear(self):
        return True

    def add_dc_stamps(self, A, z, x, iidx):
        z[self.n1] = z[self.n1] + self.I() 
        z[self.n2] = z[self.n2] - self.I()

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

