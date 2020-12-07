class IdealHarmonicFilter():
    
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.freq = float(value)

    def get_num_vsources(self):
        return 0

    def is_nonlinear(self):
        return False

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        pass

    def add_ac_stamps(self, A, z, x, iidx, freq):
        g = 1e9 if freq == self.freq else 0.
        A[self.n1][self.n1] = A[self.n1][self.n1] + g
        A[self.n2][self.n2] = A[self.n2][self.n2] + g
        A[self.n1][self.n2] = A[self.n1][self.n2] - g
        A[self.n2][self.n1] = A[self.n2][self.n1] - g

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        self.add_dc_stamps(A, z, x, iidx)

    def __str__(self):
        return 'IdealHarmonicFilter: {}\nNodes = {} -> {}\nValue = {}\n'.format(self.name, self.n1, self.n2, self.freq)
