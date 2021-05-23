import numpy as np

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
        
    def add_mthb_stamps(self, Y, S, freq, freqidx):
        if freqidx == 0:
            n1 = (self.n1 - 1) * S
            n2 = (self.n2 - 1) * S
            n3 = (self.n3 - 1) * S
            n4 = (self.n4 - 1) * S
            if (self.n1 > 0 and self.n2 > 0):
                Y[n1,n2] += 1
                Y[n2,n1] -= 1
            if (self.n1 > 0 and self.n3 > 0):
                Y[n1,n3] -= 1
                Y[n3,n1] += 1
            if (self.n2 > 0 and self.n4 > 0):
                Y[n2,n4] += 1
                Y[n4,n2] -= 1
            if (self.n3 > 0 and self.n4 > 0):
                Y[n3,n4] -= 1
                Y[n4,n3] += 1
        else:
            n1 = (self.n1 - 1) * S + 2 * (freqidx - 1) + 1
            n2 = (self.n2 - 1) * S + 2 * (freqidx - 1) + 1
            n3 = (self.n3 - 1) * S + 2 * (freqidx - 1) + 1
            n4 = (self.n4 - 1) * S + 2 * (freqidx - 1) + 1
            Ymnk = np.array([[1, 0],[0, 1]])
            if (self.n1 > 0 and self.n2 > 0):
                Y[n1:n1+2,n2:n2+2] += Ymnk
                Y[n2:n2+2,n1:n1+2] -= Ymnk
            if (self.n1 > 0 and self.n3 > 0):
                Y[n1:n1+2,n3:n3+2] -= Ymnk
                Y[n3:n3+2,n1:n1+2] += Ymnk
            if (self.n2 > 0 and self.n4 > 0):
                Y[n2:n2+2,n4:n4+2] += Ymnk
                Y[n4:n4+2,n2:n2+2] -= Ymnk
            if (self.n3 > 0 and self.n4 > 0):
                Y[n3:n3+2,n4:n4+2] -= Ymnk
                Y[n4:n4+2,n3:n3+2] += Ymnk

    def __str__(self):
        return 'Gyrator: {}\nNodes = {} -> {} and {}->{}\nG = {}\ndelay={}\n'.format(self.name, self.n1, self.n2, self.n3, self.n4, self.G)

