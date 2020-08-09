import numpy as np

import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

from yarf.Devices import *
from yarf.Analyses import DC

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('[%(levelname)s]: %(name)s: %(message)s')

file_handler = logging.FileHandler('AC.log')
file_handler.setFormatter(formatter)

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

# logger.addHandler(file_handler) # enable file logging
logger.addHandler(stream_handler) # enable console logging

options = dict()
options['is_sparse'] = False
options['max_iterations'] = 150
options['gmin'] = 1e-12

# convergence parametes
options['reltol'] = 1e-3
options['vabstol'] = 1e-6
options['iabstol'] = 1e-12

class AC():
    
    def __init__(self, name, start, stop, numpts=10, stepsize=None, sweeptype='linear'):
        self.name = name

        self.xac = [] # list of complex solutions
        self.xdc = None

        self.start = start
        self.stop = stop
        self.numpts = numpts
        self.stepsize = stepsize
        self.sweeptype = sweeptype
        
        self.n = 0 # number of nodes in the circuit with 'gnd'
        self.m = 0 # number of independent voltage sources
        self.devs = [] # list of devices in the circuit
        self.iidx = {} # map indep. vsource idx in the MNA to a device
        
        self.options = options.copy() # DC simulation options

        if sweeptype == 'linear':
            if stepsize is not None:
                self.freqs = np.arange(start, stop, step)
            else:
                self.freqs = np.linspace(start, stop, numpts)
        elif sweeptype == 'logarithm':
            self.freqs = np.geomspace(start, stop, numpts)
        else:
            logger.warning('Failed to calculate the frequencies vector!')
            self.freqs = np.array([start, stop])

    def get_dc_solution(self):
        return self.xdc

    def get_ac_solution(self):
        return self.xac

    def get_freqs(self):
        return self.freqs

    def run(self, y, x0=None):
        # Here we go!
        logger.info('Starting AC analysis.')

        # get netlist parameters and data structures
        self.n = y.get_n()
        self.m = y.get_m('ac')
        self.devs = y.get_devices()
        self.iidx = y.get_mna_extra_rows_dict('ac')

        # create MNA matrices for linear devices
        A = np.zeros((self.n+self.m, self.n+self.m), dtype=complex)
        z = np.zeros((self.n+self.m, 1), dtype=complex)

        # need to map the dc solution to the ac vector
        # because the inductor has different number of
        # rows for AC and DC simulations (or maybe not!)
        if x0 is None:
            dc = DC(self.name + 'DC')
            self.xdc = dc.run(y)
        else:
            self.xdc = x0

        self.xac.clear()
        for freq in self.freqs:
            # refresh the MNA matrices
            A[:][:] = 0
            z[:] = 0

            # populate the matrices A and z with the devices stamps
            for dev in self.devs:
                idx = self.iidx[dev] if dev in self.iidx else None
                dev.add_ac_stamps(A, z, self.xdc, idx, freq)

            # adding gmin to all nodes
            A = A + np.eye(len(A)) * self.options['gmin']

            # solve complex linear system
            xac, issolved = self.solve_ac_linear(A, z)

            # if the system is solved, add the solution to output,
            # otherwise add zeroes as solution and issue error
            if issolved:    
                self.xac.append(xac)
            else:
                self.xac.append(np.zeros((len(A)-1,1)))
                logger.error('Failed to solve AC for frequency {}!', freq)

        return self.xac

    def solve_ac_linear(self, A, z):
        if self.options['is_sparse'] == True:
            An = scipy.sparse.csc_matrix(A[1:,1:])
            lu = scipy.sparse.linalg.splu(An)
            xac = lu.solve(z[1:])
        else:
            lu, piv = scipy.linalg.lu_factor(A[1:,1:])
            xac = scipy.linalg.lu_solve((lu, piv), z[1:])

        return xac, not np.isnan(np.sum(xac))

    # TODO: remove this section, get only the convergence criteria
    def solve_dc_nonlinear(self, A, z, x0):
        # get the configuration parameters
        reltol = self.options['reltol']
        vabstol = self.options['vabstol']
        iabstol = self.options['iabstol']
        maxiter = self.options['max_iterations']

        xk = x0.copy()
        Anl = np.zeros(A.shape)
        znl = np.zeros(z.shape)
        converged = False
        k = 0
        while (not converged) and (k < maxiter):
            # refresh nonlinear MNA matrices
            Anl[:,:] = 0.0
            znl[:] = 0.0

            # add nonlinear devices to the MNA system
            for dev in self.nonlin_devs:
                idx = self.iidx[dev] if dev in self.iidx else None
                dev.add_dc_stamps(Anl, znl, xk, idx)

            # index slicing is used to remove the 'gnd' node
            # Anl and znl add the nonlinear element stamps
            An = A[1:,1:] + Anl[1:,1:]
            zn = z[1:] + znl[1:]

            # TODO: use sparse methods if 'is_sparse' is set
            #       add exception handling here (?)
            lu, piv = scipy.linalg.lu_factor(An)
            self.x  = scipy.linalg.lu_solve((lu, piv), zn)

            if np.isnan(np.sum(self.x)):
                logger.debug('Failed to resolve linear system! Solution has NaN ...') # why?
                return False

            dx = self.x - xk

            # check if voltages converged
            vconverged = True
            for i in range(0, self.n-1):
                if np.abs(dx[i,0]) > reltol * np.abs(xk[i,0]) + vabstol:
                    vconverged = False
            
            # check if currents converged
            iconverged = True
            for i in range(self.n, len(dx)):
                if np.abs(dx[i,0]) > reltol * np.abs(xk[i,0]) + iabstol:
                    iconverged = False

            if vconverged and iconverged:
                converged = True
            else:
                xk = self.x
                k = k + 1

            logger.debug('A:\n{}\nAnl:\n{}\nz:\n{}\nznl\n{}'.format(A, Anl[1:,1:], z, znl[1:]))
            logger.debug('xk:\n{}\nx:\n{}'.format(xk, self.x))

        return converged


