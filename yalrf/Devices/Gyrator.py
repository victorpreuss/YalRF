class Gyrator():
    
    def __init__(self, name, n1, n2, n3, n4, G=1):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.n3   = n3
        self.n4   = n4
        self.G    = float(G)

    def get_num_vsources(self):
        return 0

    def is_nonlinear(self):
        return False

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        A[self.n1,self.n2] = A[self.n1,self.n2] + self.G
        A[self.n1,self.n3] = A[self.n1,self.n3] - self.G
        A[self.n2,self.n1] = A[self.n2,self.n1] - self.G
        A[self.n2,self.n4] = A[self.n2,self.n4] + self.G
        A[self.n3,self.n1] = A[self.n3,self.n1] + self.G
        A[self.n3,self.n4] = A[self.n3,self.n4] - self.G
        A[self.n4,self.n2] = A[self.n4,self.n2] - self.G
        A[self.n4,self.n3] = A[self.n4,self.n3] + self.G

    def add_ac_stamps(self, A, z, x, iidx, freq):
        self.add_dc_stamps(A, z, x, iidx)

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        self.add_dc_stamps(A, z, x, iidx)

    def __str__(self):
        return 'VCCS: {}\nNodes = {} -> {} and {}->{}\nG = {}\ndelay={}\n'.format(self.name, self.n1, self.n2, self.n3, self.n4, self.G, self.tau)

