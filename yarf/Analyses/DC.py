import numpy as np

import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

from yarf.Devices import *
from yarf.Utils import dc_logger as logger

options = dict()
options['is_sparse'] = False
options['max_iterations'] = 150
options['gmin'] = 1e-12

# convergence parametes
options['reltol'] = 1e-3
options['vabstol'] = 1e-6
options['iabstol'] = 1e-12

# continuation methods
options['use_gmin_stepping'] = True
options['use_source_stepping'] = True

# parameters for the gmin stepping algorithm
options['gmin_max'] = 0.01    # initial conductance
options['gmin_min'] = 1e-12   # convergence criteria
options['gmin_rate'] = 0.05   # rate of convergence
options['gmin_maxiter'] = 50  # maximum number of iterations

# parameters for the source stepping algorithm
options['alpha_min'] = 1e-3
options['alpha_delta'] = 1e-3
options['source_maxiter'] = 50

class DC():
    
    def __init__(self, name, nodeset=None):
        self.name = name

        self.nodeset = nodeset
        self.x = None
        
        self.n = 0            # number of uniquely named nodes including 'gnd'
        self.m = 0            # number of independent voltage sources
        self.iidx = {}        # maps an indep. vsource idx in the MNA to a device
        self.lin_devs = []    # list of linear devices
        self.nonlin_devs = [] # list of nonlinear devices
        
        self.options = options.copy() # DC simulation options

    def get_dc_solution(self):
        return self.x

    def run(self, y, x0=None):
        # get netlist parameters and data structures
        self.n = y.get_n()
        self.m = y.get_m('dc')
        self.lin_devs = y.get_linear_devices()
        self.nonlin_devs = y.get_nonlinear_devices()
        self.iidx = y.get_mna_extra_rows_dict('dc')

        # create MNA matrices for linear devices
        A = np.zeros((self.n+self.m, self.n+self.m))
        z = np.zeros((self.n+self.m, 1))

        # TODO: create initial guess vector from nodeset (forget x0)
        x0 = np.zeros((len(A)-1, 1)) if x0 is None else x0

        # some devices have initialization routines
        for dev in y.get_devices():
            dev.init()

        # Here we go!
        logger.info('Starting DC analysis.')

        # populate the matrices A and z with the linear devices stamps
        for dev in self.lin_devs:
            idx = self.iidx[dev] if dev in self.iidx else None
            dev.add_dc_stamps(A, z, None, idx)

        # TODO: add gmin only at problematic nodes such as middle of two
        #       capacitors or in parallel to pn junctions (high conductance) to
        #       prevent the occurence of a singular matrix.
        A = A + np.eye(len(A)) * self.options['gmin']

        # if there is not a nonlinear device, simply solve the linear system
        if not self.nonlin_devs:
            logger.info('Starting linear DC solver ...')
            self.x, issolved = self.solve_linear(A[1:,1:], z[1:])
            if issolved:
                logger.info('Finished DC analysis.')
                return self.x

        # solve nonlinear DC analysis using Newton-Raphson
        logger.info('Starting nonlinear DC solver ...')
        issolved = self.solve_dc_nonlinear(A, z, x0)
        if issolved:
            logger.info('Finished DC analysis.')
            return self.x

        # solve DC analysis using the gmin stepping continuation method
        if self.options['use_gmin_stepping'] == True:
            logger.info('Previous solver failed! Using gmin stepping ...')
            issolved = self.solve_dc_nonlinear_using_gmin_stepping(A, z, x0)
            if issolved:
                logger.info('Finished DC analysis.')
                return self.x

        # solve DC analysis using the source stepping continuation method
        if self.options['use_source_stepping'] == True:
            logger.info('Previous solver failed! Using source stepping ...')
            issolved = self.solve_dc_nonlinear_using_source_stepping(A, z, x0)
            if issolved:
                logger.info('Finished DC analysis.')
                return self.x

        # TODO: implement more continuation/homotopy techniques here

        logger.error('DC analysis failed to converge!')
        return None
    
    def solve_linear(self, A, z, is_sparse=False):
        if is_sparse == True:
            An = scipy.sparse.csc_matrix(A)
            lu = scipy.sparse.linalg.splu(An)
            x = lu.solve(z)
        else:
            lu, piv = scipy.linalg.lu_factor(A)
            x = scipy.linalg.lu_solve((lu, piv), z)

        return x, not np.isnan(np.sum(x))

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
            # refresh matrices
            Anl[:,:] = 0.0
            znl[:] = 0.0

            # add nonlinear element stamps
            for dev in self.nonlin_devs:
                idx = self.iidx[dev] if dev in self.iidx else None
                dev.add_dc_stamps(Anl, znl, xk, idx)

            # index slicing is used to remove the 'gnd' node
            # An is the Jacobian matrix of the Newton-Raphson iteration
            An = A[1:,1:] + Anl[1:,1:]
            zn = z[1:] + znl[1:]

            # solve linear system
            self.x, issolved = self.solve_linear(An, zn, self.options['is_sparse'])

            if not issolved:
                logger.debug('Failed to resolve linear system! Solution has NaN ...')
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

            logger.debug('A:\n{}\nAnl:\n{}\nz:\n{}\nznl\n{}'.format(A, Anl, z, znl))
            logger.debug('xk:\n{}\nx:\n{}'.format(xk, self.x))

        return converged

    def solve_dc_nonlinear_using_gmin_stepping(self, A, z, x0):
        # get the configuration parameters
        gmin_max = self.options['gmin_max']
        gmin_min = self.options['gmin_min']
        gmin_rate = self.options['gmin_rate']
        gmin_maxiter = self.options['gmin_maxiter']

        gmin = gmin_max
        xprev = x0.copy()
        self.x = x0.copy()
        converged = False
        k = 0
        while (not converged) and (k < gmin_maxiter):
            # add (decreasing) conductance from every node to 'gnd'
            An = A + np.eye(len(A)) * gmin

            # solve for this 'augmented' circuit using previous solution
            issolved = self.solve_dc_nonlinear(An, z, self.x)

            # check if convergence is reached
            if issolved:
                xprev = self.x # backup of good solution

                if gmin < gmin_min:
                    converged = True
                else:
                    gmin = gmin * gmin_rate
                    gmin_rate = gmin_rate / 2
            else:
                self.x = xprev # restore backup of good solution

                # decrease rate of gmin stepping if convergence fails
                gmin = gmin / gmin_rate
                gmin_rate = 2 * gmin_rate

                if gmin > gmin_max:
                    converged = False
                    break
            k = k + 1

        return converged

    def solve_dc_nonlinear_using_source_stepping(self, A, z, x0):
        # get the configuration parameters
        alpha_min = self.options['alpha_min']
        alpha_delta = self.options['alpha_delta']
        maxiter = self.options['source_maxiter']

        alpha = alpha_min
        xprev = x0.copy()
        self.x = x0.copy()
        converged = False
        k = 0
        while (not converged) and (k < maxiter):
            # sequentially increase the voltage and current supplies
            zn = alpha * z

            # solve for this version of the circuit using previous solution
            issolved = self.solve_dc_nonlinear(A, zn, self.x)

            # check if convergence is reached
            if issolved:
                # backup of good solution
                xprev = self.x 

                if alpha >= 1.0:
                    converged = True
                else:
                    alpha = alpha + alpha_delta
                    alpha_delta = 2 * alpha_delta
                    if alpha >= 1.0:
                        alpha = 1.0
            else:
                # restore backup of good solution
                self.x = xprev 

                # decrease rate of source stepping if convergence fails
                alpha = alpha - alpha_delta / 2
                alpha_delta = alpha_delta / 2

                if alpha > alpha_min:
                    converged = False
                    break
            k = k + 1

        return converged

