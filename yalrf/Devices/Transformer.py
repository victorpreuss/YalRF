import numpy as np

class Transformer():
    
    def __init__(self, name, n1, n2, n3, n4, T=1):
        self.name = name
        self.n1   = n1
        self.n2   = n2
        self.n3   = n3
        self.n4   = n4
        self.T    = float(T)

    def get_num_vsources(self):
        return 1

    def is_nonlinear(self):
        return False

    def init(self):
        pass

    def add_dc_stamps(self, A, z, x, iidx):
        A[self.n1,iidx] = A[self.n1,iidx] - 1.
        A[self.n2,iidx] = A[self.n2,iidx] + self.T
        A[self.n3,iidx] = A[self.n3,iidx] - self.T
        A[self.n4,iidx] = A[self.n4,iidx] + 1.
        A[iidx,self.n1] = A[iidx,self.n1] + 1.
        A[iidx,self.n2] = A[iidx,self.n2] - self.T
        A[iidx,self.n3] = A[iidx,self.n3] + self.T
        A[iidx,self.n4] = A[iidx,self.n4] - 1.

    def add_ac_stamps(self, A, z, x, iidx, freq):
        self.add_dc_stamps(A, z, x, iidx)

    def add_tran_stamps(self, A, z, x, iidx, xt, t, tstep):
        self.add_dc_stamps(A, z, x, iidx)
        
    def add_mthb_stamps(self, Y, S, freq, freqidx):
        # TODO
        pass

    def __str__(self):
        return 'Transformer: {}\nNodes = {} -> {} and {}->{}\nT = {}\n'.format(self.name, self.n1, self.n2, self.n3, self.n4, self.T)

