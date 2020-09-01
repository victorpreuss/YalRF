class VoltageControlledVoltageSource():
    
    def __init__(self, name, n1, n2, n3, n4, G=1, tau=None):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.n3   = n3
        self.n4   = n4
        self.G    = G
        self.tau  = tau

    def get_num_vsources(self, analysis):
        return 1

    def is_nonlinear(self):
        return False

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        A[iidx][self.n1] = +self.G
        A[iidx][self.n2] = -1.0 
        A[iidx][self.n3] = +1.0 
        A[iidx][self.n4] = -self.G
        A[self.n2][iidx] = -1.0
        A[self.n3][iidx] = +1.0

    def __str__(self):
        return 'VCVS: {}\nNodes = {} -> {} and {}->{}\nG = {}\ndelay={}\n'.format(self.name, self.n1, self.n2, self.n3, self.n4, self.G, self.tau)

