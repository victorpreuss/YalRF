from .Devices import *
from .Analyses import *
from .Utils import yalrf_logger as logger
from .Netlist import Netlist

import logging
logger.setLevel(logging.WARNING)

class YalRF(Netlist):
    """
    Class containing the YalRF API for circuit simulation.
    
    This class is initialized with only the circuit name. It contains methods to
    create and manipulate netlists, add analyses to be performed in the circuit,
    and manipulate the solution data from analyses.

    """
    def __init__(self, name):
        """
        Class constructor.

        Creates an empty netlist with only the 'gnd' node.

        Parameters
        ----------
        name : str
            Name of the circuit.

        """
        self.name = name

        # initialize circuit netlist
        Netlist.__init__(self, name)

        # list of analysis objects (DC, AC, Transient, ...)
        self.analyses = {}

    def add_tran_analysis(self, name, tstop, maxtstep=None, tstart=0, uic=False):
        """
        Create and add a Transient analysis.

        Parameters
        ----------
        name : str
            Name for the analysis object.
        tstop : float
            Stop time of the simulation.
        maxtstep : float
            Maximum timestep value.
        tstart : float
            Time to start saving simulation output data.
        uic : bool
            Use provided initial conditions for capacitors and inductors
            instead of calculating the DC bias point at the first step.

        Returns
        -------
        :class:`Transient`
            Reference to the created Transient object. This allows the user to change
            internal parameters of the instance before running it.

        """
        if name not in self.analyses:
            tran = Transient(name, tstop, maxtstep, tstart, uic)
            self.analyses[name] = tran
            return tran
        else:
            logger.warning('Analysis name \'{}\' already taken!'.format(name))
            return None

    def add_ac_analysis(self, name, start, stop, numpts=10, stepsize=None, sweeptype='linear'):
        """
        Create and add an AC analysis.

        Parameters
        ----------
        name : str
            Name for the analysis object.
        start : float
            Start frequency.
        stop : float
            Stop frequency.
        numpts : int
            Number of points in the sweep between start and stop.
        stepsize : float
            Difference between two subsequent frequency points. Only available
            in linear sweep. It defines the number of points in the sweep.
        sweeptype : str
            Type of sweep to be performed. Possible values are: linear, logarithm.

        Returns
        -------
        :class:`AC`
            Reference to the created AC object. This allows the user to change
            internal parameters of the instance before running it.

        """
        if name not in self.analyses:
            ac = AC(name, start, stop, numpts, stepsize, sweeptype)
            self.analyses[name] = ac
            return ac
        else:
            logger.warning('Analysis name \'{}\' already taken!'.format(name))
            return None

    def add_dc_analysis(self, name):
        """
        Create and add a DC analysis with an optional nodeset.

        Parameters
        ----------
        name : str
            Name for the analysis object.

        Returns
        -------
        :class:`DC`
            Reference to the created DC object. This allows the user to change
            internal parameters of the instance before running it.

        """
        if name not in self.analyses:
            dc = DC(name)
            self.analyses[name] = dc
            return dc
        else:
            logger.warning('Analysis name \'{}\' already taken!'.format(name))
            return None

    def run(self, name, x0=None, nodeset=None):
        """
        Run analysis with the requested name using initial condition.

        Parameters
        ----------
        name : str
            Name for the analysis object to run.
        x0 : :class:`numpy.ndarray`
            DC solution to be used as initial condition or operating point for
            the analysis. The DC analysis use this as an initial guess.
        nodeset : dict
            Dictionary with node names and guesses for their initial condition.

        Returns
        -------
        list or :class:`numpy.ndarray`
            Solution data. Depending on the analysis (AC, DC, Transient, ...)
            it might be stored in different data structures.

        """
        if name in self.analyses:
            sol = self.analyses[name].run(self, x0, nodeset)
            return sol
        else:
            logger.warning('Unknown analysis name: {}!'.format(name))
            return None

    def get_analysis(self, name):
        """
        Return reference to an analysis object with the requested name.

        Parameters
        ----------
        name : str
            Name of the analysis to get instance.

        Returns
        -------
        :class:`Analyses`
            Analysis with the requested name.

        """
        if name in self.analyses:
            a = self.analyses[name]
            return a
        else:
            logger.warning('Unknown analysis name: {}!'.format(name))
            return None

    def get_voltage(self, analysis, node):
        """
        Return the voltage result for a determined node and analysis.

        Parameters
        ----------
        analysis : str
            Name of the analysis to get the results from.
        node : str
            Name of the node to consult the result voltage.

        Returns
        -------
        :class:`numpy.ndarray` or float
            Voltage result for a determined node and analysis.

        """
        a = self.get_analysis(analysis)
        if isinstance(a, DC):
            v = a.get_dc_solution()[self.get_voltage_idx(node), 0]
            return v
        elif isinstance(a, AC):
            v = a.get_ac_solution()[:, self.get_voltage_idx(node)]
            return v
        elif isinstance(a, Transient):
            v = a.get_tran_solution()[:, self.get_voltage_idx(node)]
            return v
        else:
            logger.warning('Unknown analysis type!')
            return None

    def get_idc(self, analysis, device):
        a = self.get_analysis(analysis)
        dev = self.get_device(device)
        if dev in a.iidx:
            i = a.get_dc_solution()[a.iidx[dev]-1, 0]
            return i
        else:
            x = a.get_dc_solution()
            i = dev.get_idc(x)
            return i

    def get_itran(self, analysis, device):
        a = self.get_analysis(analysis)
        dev = self.get_device(device)
        if isinstance(a, Transient):
            if dev in a.iidx:
                i = a.get_tran_solution()[:, a.iidx[dev]-1, 0]
                return i
            else:
                x = a.get_tran_solution()
                i = dev.get_itran(x)
                return i
        else:
            logger.warning('Not a transient simulation!')
            return None

    def get_time(self, analysis):
        """
        Return the simulated time array of a Transient analysis.

        Parameters
        ----------
        analysis : str
            Name of the analysis to get the time array from.
        
        Returns
        -------
        :class:`numpy.ndarray`
            Time array.

        """
        a = self.get_analysis(analysis)
        if isinstance(a, Transient):
            return a.get_time()
        else:
            logger.warning('Analysis doesn\'t have a time array!')
            return None

    def get_freqs(self, analysis):
        """
        Return the frequency array of an AC analysis.

        Parameters
        ----------
        analysis : str
            Name of the analysis to get the frequencies from.
        
        Returns
        -------
        :class:`numpy.ndarray`
            Frequency array.

        """
        a = self.get_analysis(analysis)
        if isinstance(a, AC):
            return a.get_freqs()
        else:
            logger.warning('Analysis doesn\'t have a frequency array!')
            return None

    def print_dc_voltages(self, name):
        """
        Print the DC voltages in the circuit.

        Parameters
        ----------
        name : str
            Name of the analysis to print the results.

        """
        if name in self.analyses:
            x = self.analyses[name].get_dc_solution()
            print('DC Voltages:')
            for i in range(1, self.get_n()):
                print('{}:\t{:0.4f} V'.format(self.node_idx_to_name[i], x[i-1,0]))
            print()
        else:
            logger.warning('Unknown analysis name: {}!'.format(name))
    
    def print_dc_currents(self, name):
        """
        Print the DC currents in the circuit.

        Parameters
        ----------
        name : str
            Name of the analysis to print the results.

        """
        if name in self.analyses:
            x = self.analyses[name].get_dc_solution()
            n = self.get_n()
            k = 0
            print('DC Currents:')
            for dev in self.devices:
                if dev.get_num_vsources() == 1:
                    print('{}:\t{:0.6f} A'.format(dev.name, x[n+k-1,0]))
                elif dev.get_num_vsources() == 2:
                    print('{}:\t{:0.6f} A (I1)'.format(dev.name, x[n+k-1,0]))
                    print('{}:\t{:0.6f} A (I2)'.format(dev.name, x[n+k,0]))
                k = k + dev.get_num_vsources()
            print()
        else:
            logger.warning('Unknown analysis name: {}!'.format(name))


    




