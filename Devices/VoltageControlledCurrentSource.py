class VoltageControlledCurrentSource():
    
    def __init__(self, name, n1, n2, n3, n4, G=1, tau=None):
        self.name = name
        self.n1    = n1
        self.n2    = n2
        self.n3    = n3
        self.n4    = n4
        self.G     = G
        self.tau   = tau

    def get_vsource(self):
        return 0

    def is_linear(self):
        return True

    def add_dc_stamps(self, A, z, x, iidx):
        A[self.n2][self.n1] = A[self.n3][self.n1] + self.G
        A[self.n2][self.n4] = A[self.n3][self.n1] - self.G
        A[self.n3][self.n1] = A[self.n3][self.n1] - self.G
        A[self.n3][self.n4] = A[self.n3][self.n1] + self.G

    def __str__(self):
        return 'VCCS: {}\nNodes = {} -> {} and {}->{}\nG = {}\ndelay={}\n'.format(self.name, self.n1, self.n2, self.n3, self.n4, self.G, self.tau)

