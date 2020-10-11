class CurrentControlledCurrentSource():
    
    def __init__(self, name, n1, n2, n3, n4, G=1, tau=0):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.n3   = n3
        self.n4   = n4
        self.G    = float(G)
        self.tau  = float(tau)

    def get_num_vsources(self):
        return 1

    def is_nonlinear(self):
        return False

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        A[iidx][self.n1] = +1.0
        A[iidx][self.n4] = -1.0
        A[self.n1][iidx] = +1.0
        A[self.n2][iidx] = +self.G
        A[self.n3][iidx] = -self.G
        A[self.n4][iidx] = -1.0

    def add_ac_stamps(self, A, z, x, iidx, freq):
        G = self.G * np.exp(-1j * 2. * np.pi * freq * self.tau)
        A[iidx][self.n1] = +1.0
        A[iidx][self.n4] = -1.0
        A[self.n1][iidx] = +1.0
        A[self.n2][iidx] = +G
        A[self.n3][iidx] = -G
        A[self.n4][iidx] = -1.0
    
    # TODO: implement transient time delay
    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        self.add_dc_stamps(A, z, x, iidx)

    def __str__(self):
        return 'CCCS: {}\nNodes = {} -> {} and {}->{}\nG = {}\ndelay={}\n'.format(self.name, self.n1, self.n2, self.n3, self.n4, self.G, self.tau)

