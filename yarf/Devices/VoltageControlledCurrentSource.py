class VoltageControlledCurrentSource():
    
    def __init__(self, name, n1, n2, n3, n4, G=1, tau=0):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.n3   = n3
        self.n4   = n4
        self.G    = float(G)
        self.tau  = float(tau)

    def get_num_vsources(self, analysis):
        return 0

    def is_nonlinear(self):
        return False

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        A[self.n2][self.n1] = A[self.n3][self.n1] + self.G
        A[self.n2][self.n4] = A[self.n3][self.n1] - self.G
        A[self.n3][self.n1] = A[self.n3][self.n1] - self.G
        A[self.n3][self.n4] = A[self.n3][self.n1] + self.G

    def add_ac_stamps(self, A, z, x, iidx, freq):
        G = self.G * np.exp(-1j * 2. * np.pi * freq * self.tau)
        A[self.n2][self.n1] = A[self.n3][self.n1] + G
        A[self.n2][self.n4] = A[self.n3][self.n1] - G
        A[self.n3][self.n1] = A[self.n3][self.n1] - G
        A[self.n3][self.n4] = A[self.n3][self.n1] + G

    def __str__(self):
        return 'VCCS: {}\nNodes = {} -> {} and {}->{}\nG = {}\ndelay={}\n'.format(self.name, self.n1, self.n2, self.n3, self.n4, self.G, self.tau)

