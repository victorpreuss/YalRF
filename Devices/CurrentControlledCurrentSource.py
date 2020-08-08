class CurrentControlledCurrentSource():
    
    def __init__(self, name, n1, n2, n3, n4, G=None, tau=None):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.n3   = n3
        self.n4   = n4
        self.G    = G
        self.tau  = tau

    def get_num_vsources(self, analysis):
        return 1

    def is_linear(self):
        return True

    def add_dc_stamps(self, A, z, x, iidx):
        A[iidx][self.n1] = +1.0
        A[iidx][self.n4] = -1.0 
        A[self.n1][iidx] = +1.0
        A[self.n2][iidx] = +self.G
        A[self.n3][iidx] = -self.G
        A[self.n4][iidx] = -1.0

    def __str__(self):
        return 'CCCS: {}\nNodes = {} -> {} and {}->{}\nG = {}\ndelay={}\n'.format(self.name, self.n1, self.n2, self.n3, self.n4, self.G, self.tau)

