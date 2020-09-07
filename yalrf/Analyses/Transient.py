import sys
import numpy as np

from yalrf.Devices import *
from yalrf.Analyses import DC
from yalrf.Analyses.Solver import solve_linear
from yalrf.Utils import tr_logger as logger

import logging
logger.setLevel(logging.INFO)

options = dict()
options['is_sparse'] = False
options['max_iterations'] = 150
options['gmin'] = 1e-12

# convergence parameters
options['reltol'] = 1e-3
options['vabstol'] = 1e-6
options['iabstol'] = 1e-12

# transient parameters
options['mintstep'] = 1e-16

class Transient():

    def __init__(self, name, tstart, tstop, maxtstep=None, uic=False):
        self.name = name

        # output data
        self.xtran = None
        self.xdc = None
        self.time = None

        # analysis parameters
        self.tstart = tstart
        self.tstop = tstop
        self.maxtstep = maxtstep
        self.uic = uic

        self.options = options.copy() # Transient simulation options

    def get_dc_solution(self):
        return self.xdc

    def get_tran_solution(self):
        return self.xtran

    def get_time(self):
        return self.time

    def run(self, y, x0=None, nodeset=None):
        # get netlist parameters and data structures
        self.n = y.get_n()
        self.m = y.get_m('tran')
        self.devs = y.get_devices()
        self.iidx = y.get_mna_extra_rows_dict('tran')

        # create MNA matrices
        A = np.zeros((self.n+self.m, self.n+self.m)) 
        z = np.zeros((self.n+self.m, 1))

        # perform DC simulation if no initial point is provided
        if x0 is None:
            dc = DC(self.name + '.DC')
            self.xdc = dc.run(y, nodeset=nodeset)
        else:
            self.xdc = x0

         # initialize devices and calculate the operating point of nonlinear elements
        for dev in self.devs:
            dev.init()
            if dev.is_nonlinear():
                dev.calc_oppoint(self.xdc)

        # Here we go!
        logger.info('Starting Transient analysis.')

        # get the configuration parameters
        reltol = self.options['reltol']
        vabstol = self.options['vabstol']
        iabstol = self.options['iabstol']
        maxiter = self.options['max_iterations']
        mintstep = self.options['mintstep']

        j = 0               # iterator
        t = 0.              # time variable
        tstep = 1e-15       # time step
        xtran = [self.xdc]  # output data
        time  = [t]         # output time array
        while t < self.tstop:
            # increment time step
            t = t + tstep

            # use last transient point as initial condition for solving the next
            xk = xtran[-1]

            converged = False
            k = 0
            while (not converged) and (k < maxiter):
                # refresh matrices
                A[:,:] = 0.0
                z[:] = 0.0

                # add transient stamps to MNA
                for dev in self.devs:
                    idx = self.iidx[dev] if dev in self.iidx else None
                    dev.add_tran_stamps(A, z, xtran, idx, t, tstep)

                # solve linear system
                x, issolved = solve_linear(A[1:,1:], z[1:])

                if not issolved:
                    logger.debug('Failed to resolve linear system! Solution has NaN ...')
                    break

                dx = x - xk

                # check if voltages converged
                vconverged = True
                for i in range(0, self.n-1):
                    if np.abs(dx[i,0]) > reltol * np.abs(xk[i,0]) + vabstol:
                        vconverged = False
                
                # check if currents converged
                iconverged = True
                for i in range(self.n-1, len(dx)):
                    if np.abs(dx[i,0]) > reltol * np.abs(xk[i,0]) + iabstol:
                        iconverged = False

                logger.debug('\nA:\n{}\nz:\n{}\n'.format(A, z))
                logger.debug('\nx:\n{}'.format(x[-1]))

                # finish algorithm if simulation converged
                if vconverged and iconverged and k >= 1:
                    converged = True
                else:
                    xk = x
                    k = k + 1

            logger.debug('The solver took {} iterations.'.format(k+1))
            logger.debug('j = ' + str(j) + '\ntstep = ' + str(tstep))

            if converged:
                # save solution
                time.append(t)
                xtran.append(x)
                j = j + 1

                # save capacitor currents
                for dev in self.devs:
                    if isinstance(dev, Capacitor):
                        dev.store_current(xtran, tstep)

                # recalculate time step
                if k < 5:
                    tstep = min(tstep * 2., self.maxtstep)
                elif k > 10:
                    tstep = tstep / 2.

            else:
                # reduce time step if NR failed to converge
                t = t - tstep
                tstep = tstep / 10.

                if tstep < mintstep:
                    logger.error('Timestep: {} is below the minimum allowed!'.format(tstep))
                    sys.exit(1)

        # transform outputs into numpy array
        self.xtran = np.stack(xtran, axis=0)
        self.time = np.array(time)

        logger.info('Finished Transient analysis.')
        return self.xtran
