from .Devices import *
from .Analyses import *
from .Utils import yalrf_logger as logger

import logging
logger.setLevel(logging.WARNING)

class YalRF():
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

        # netlist related attributes
        self.devices = []                   # list of all the devices in the netlist
        self.node_name_to_idx = {'gnd': 0}  # dictionary to associate a node name (string) to its index in the netlist
        self.node_idx_to_name = ['gnd']     # list to associate a node index in the netlist to its node name (string)

        # analysis related attributes
        self.analyses = {}                  # list of analysis objects (DC, AC, Transient, ...)

    def add_tran_analysis(self, name, tstart, tstop, maxtstep=None, uic=False):
        """
        Create and add a Transient analysis.

        Parameters
        ----------
        name : str
            Name for the analysis object.
        tstart : float
            Time to start saving simulation output data.
        tstop : float
            Stop time of the simulation.
        maxtstep : float
            Maximum timestep value.
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
            tran = Transient(name, tstart, tstop, maxtstep, uic)
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

    def add_resistor(self, name, n1, n2, value):
        """
        Add resistor to the netlist.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the resistor.
        n2 : str
            Negative node of the resistor.
        value : float
            Resistance value in ohms.

        Returns
        -------
        :class:`Resistor`
            Reference to the created Resistor object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        res = Resistor(name, n1, n2, value)
        self.devices.append(res)
        return res
        
    def add_capacitor(self, name, n1, n2, value):
        """
        Add capacitor to the netlist.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the capacitor.
        n2 : str
            Negative node of the capacitor.
        value : float
            Capacitance value in Farads.

        Returns
        -------
        :class:`Capacitor`
            Reference to the created Capacitor object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        cap = Capacitor(name, n1, n2, value)
        self.devices.append(cap)
        return cap
        
    def add_inductor(self, name, n1, n2, value):
        """
        Add inductor to the netlist.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the inductor.
        n2 : str
            Negative node of the inductor.
        value : float
            Inductance value in Henrys.

        Returns
        -------
        :class:`Inductor`
            Reference to the created Inductor object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        ind = Inductor(name, n1, n2, value)
        self.devices.append(ind)
        return ind

    def add_vdc(self, name, n1, n2, dc):
        """
        Add DC voltage source to the netlist.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the voltage source.
        n2 : str
            Negative node of the voltage source.
        dc : float
            DC voltage value in Volts.

        Returns
        -------
        :class:`VoltageSource`
            Reference to the created VoltageSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        vdc = VoltageSource(name, n1, n2, vtype='dc', dc=dc)
        self.devices.append(vdc)
        return vdc

    def add_idc(self, name, n1, n2, dc):
        """
        Add DC current source to the netlist.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the current source.
        n2 : str
            Negative node of the current source.
        dc : float
            DC current value in Amperes.

        Returns
        -------
        :class:`CurrentSource`
            Reference to the created CurrentSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        idc = CurrentSource(name, n1, n2, itype='dc', dc=dc)
        self.devices.append(idc)
        return idc

    def add_vac(self, name, n1, n2, ac, phase=0):
        """
        Add AC voltage source to the netlist.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the voltage source.
        n2 : str
            Negative node of the voltage source.
        ac : float
            AC voltage value in Volts.
        phase : float
            Phase of the AC voltage source

        Returns
        -------
        :class:`VoltageSource`
            Reference to the created VoltageSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        vac = VoltageSource(name, n1, n2, vtype='ac', ac=ac, phase=phase)
        self.devices.append(vac)
        return vac

    def add_iac(self, name, n1, n2, ac, phase=0):
        """
        Add AC current source to the netlist.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the current source.
        n2 : str
            Negative node of the current source.
        ac : float
            AC current value in Amperes.
        phase : float
            Phase of the AC current source

        Returns
        -------
        :class:`CurrentSource`
            Reference to the created CurrentSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        iac = CurrentSource(name, n1, n2, itype='ac', ac=ac, phase=phase)
        self.devices.append(iac)
        return iac

    def add_vsource(self, name, n1, n2, dc, ac, phase=0):
        """
        Add a voltage source to the netlist with DC and AC values.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the voltage source.
        n2 : str
            Negative node of the voltage source.
        dc : float
            DC voltage value in Volts.
        ac : float
            AC voltage value in Volts.
        phase : float
            Phase of the AC voltage source

        Returns
        -------
        :class:`VoltageSource`
            Reference to the created VoltageSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        vsource = VoltageSource(name, n1, n2, dc, ac, phase, vtype='both')
        self.devices.append(vsource)
        return vsource

    def add_isource(self, name, n1, n2, dc, ac, phase=0):
        """
        Add a current source to the netlist with DC and AC values.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the current source.
        n2 : str
            Negative node of the current source.
        dc : float
            DC current value in Amperes.
        ac : float
            AC current value in Amperes.
        phase : float
            Phase of the AC current source

        Returns
        -------
        :class:`CurrentSource`
            Reference to the created CurrentSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        isource = CurrentSource(name, n1, n2, dc, ac, phase, itype='both')
        self.devices.append(isource)
        return isource

    def add_vstep(self, name, n1, n2, dc, tstart, tstop=1e6):
        """
        Add a voltage step source to the netlist. This source is used to apply
        voltage steps at Transient simulations. It is 0V fot DC and AC analysis.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the voltage source.
        n2 : str
            Negative node of the voltage source.
        dc : float
            DC voltage value in Volts.
        tstart : float
            Time at which the voltage step is applied.
        tstop : float
            Time at which the voltage step ends.

        Returns
        -------
        :class:`TransientVoltageSource`
            Reference to the created VoltageSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        vstep = TransientVoltageSource(name, n1, n2, dc, tstart, tstop, vtype='step')
        self.devices.append(vstep)
        return vstep

    def add_vcvs(self, name, n1, n2, n3, n4, G, tau=0):
        """
        Add a voltage controlled voltage source to the netlist.

        n1 o----- +   + -->--o n2
                 Vin Vout
        n4 o----- -   - -----o n3

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Node 1.
        n2 : str
            Node 2.
        n3 : str
            Node 3.
        n4 : str
            Node 4.
        G : float
            Gain/Transfer factor.
        tau : float
            Delay time.

        Returns
        -------
        :class:`VoltageControlledVoltageSource`
            Reference to the created VoltageControlledVoltageSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        n4 = self.add_node(n4)

        vcvs = VoltageControlledVoltageSource(name, n1, n2, n3, n4, G, tau)
        self.devices.append(vcvs)
        return vcvs

    def add_vccs(self, name, n1, n2, n3, n4, G, tau=0):
        """
        Add a voltage controlled current source to the netlist.

        n1 o----- +   ---<--o n2
                 Vin  | Iout
        n4 o----- -   --->--o n3

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Node 1.
        n2 : str
            Node 2.
        n3 : str
            Node 3.
        n4 : str
            Node 4.
        G : float
            Gain/Transfer factor.
        tau : float
            Delay time.

        Returns
        -------
        :class:`VoltageControlledCurrentSource`
            Reference to the created VoltageControlledCurrentSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        n4 = self.add_node(n4)

        vccs = VoltageControlledCurrentSource(name, n1, n2, n3, n4, G, tau)
        self.devices.append(vccs)
        return vccs
        
    def add_ccvs(self, name, n1, n2, n3, n4, G, tau=0):
        """
        Add a current controlled voltage source to the netlist.

        n1 o--->--   + -->--o n2
             Iin |  Vout
        n4 o---<--   - -----o n3

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Node 1.
        n2 : str
            Node 2.
        n3 : str
            Node 3.
        n4 : str
            Node 4.
        G : float
            Gain/Transfer factor.
        tau : float
            Delay time.

        Returns
        -------
        :class:`CurrentControlledVoltageSource`
            Reference to the created CurrentControlledVoltageSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        n4 = self.add_node(n4)

        ccvs = CurrentControlledVoltageSource(name, n1, n2, n3, n4, G, tau)
        self.devices.append(ccvs)
        return ccvs
        
    def add_cccs(self, name, n1, n2, n3, n4, G, tau=0):
        """
        Add a current controlled current source to the netlist.

        n1 o--->--   --<---o n2
             Iin |   | Iout
        n4 o---<--   -->---o n3

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Node 1.
        n2 : str
            Node 2.
        n3 : str
            Node 3.
        n4 : str
            Node 4.
        G : float
            Gain/Transfer factor.
        tau : float
            Delay time.

        Returns
        -------
        :class:`CurrentControlledCurrentSource`
            Reference to the created CurrentControlledCurrentSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        n4 = self.add_node(n4)

        cccs = CurrentControlledCurrentSource(name, n1, n2, n3, n4, G, tau)
        self.devices.append(cccs)
        return cccs

    def add_diode(self, name, n1, n2):
        """
        Add a diode to the netlist.

        The instance of the added diode is returned for the user. A dictionary
        containing all the parameters for the diode model can be accessed and
        modified before the simulation is ran, to configure the model. Access
        the file Devices/Diode.py for a list of available parameters.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Anode (+).
        n2 : str
            Cathode (-).

        Returns
        -------
        :class:`Diode`
            Reference to the created Diode object. Useful to edit the model options.

        Examples
        --------

        >>> d1 = y.add_diode('D1', 'n1', 'n2')
        >>> d1.options['Is'] = 1e-15

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        diode = Diode(name, n1, n2)
        self.devices.append(diode)
        return diode

    def add_bjt(self, name, n1, n2, n3, n4='gnd'):
        """
        Add a BJT to the netlist.

        The instance of the added bjt is returned for the user. A dictionary
        containing all the parameters for the bjt model can be accessed and
        modified before the simulation is ran, to configure the model. Access
        the file Devices/BJT.py for the a list of available parameters.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Base (B).
        n2 : str
            Colector (C).
        n3: str
            Emitter(E).
        n4: str
            Substrate (S).

        Returns
        -------
        :class:`BJT`
            Reference to the created BJT object. Useful to edit the model options.

        Examples
        --------

        >>> q1 = y.add_bjt('Q1', 'n1', 'n2', 'n3')
        >>> q1.options['Bf'] = 500

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        n4 = self.add_node(n4)
        
        bjt = BJT(name, n1, n2, n3)
        self.devices.append(bjt)
        return bjt

    # TODO:
    # def add_mosfet(self, name, n1, n2, n3, n4='gnd'):
    # def add_iprobe(self, name, n1, n2):
    # def add_opamp(self, name, n1, n2, n3):
    # def add_transformer(self, name, n1, n2, n3, n4, T):
    # def add_subcircuit(self, name, YalRF object):
    # def remove_device(self, name):
    # def remove_analysis(self, name):

    def add_node(self, n):
        """
        Add node and returns its index in the netlist.

        This method checks if a node name (str) already exists in the netlist.
        If it does, returns the previous index assigned for this node. Otherwise
        it adds the new node name to the netlist and assign an index to it.

        Parameters
        ----------
        n : str
            Name of the node to have an index assigned.

        Returns
        -------
        int
            Index of the node in the netlist.

        """
        if n not in self.node_name_to_idx:
            self.node_name_to_idx[n] = len(self.node_name_to_idx)
            self.node_idx_to_name.append(n)
        return self.node_name_to_idx[n]

    def is_nonlinear(self):
        """Return True if netlist has a nonlinear device."""
        for dev in self.devices:
            if dev.is_nonlinear():
                return True
        return True

    def get_n(self):
        """Return number of uniquely named nodes in the netlist."""
        return len(self.node_idx_to_name)

    def get_m(self, analysis):
        """
        Return number of independent voltage sources for an analysis.

        The number of independent voltage sources needed are important for
        the MNA algorithm, since it determines the number of additional rows
        required in the matrix system. This number depends on the companion
        models used by each device for a determined analysis.

        Parameters
        ----------
        analysis : str
            Name of the analysis.

        Returns
        -------
        int
            Number of independent voltage sources in an analysis.

        """
        m = 0
        for dev in self.devices:
            m = m + dev.get_num_vsources(analysis)
        return m

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

    def get_device(self, name):
        """
        Return reference to a device object with the requested name.

        Parameters
        ----------
        name : str
            Name of the device to get instance.

        Returns
        -------
        :class:`Devices`
            Device with the requested name.

        """
        for dev in self.devices:
            if dev.name == name:
                return dev
        logger.warning('Unknown device name: {}!'.format(name))
        return None

    def get_devices(self):
        """Return list of devices in the netlist."""
        return self.devices

    def get_linear_devices(self):
        """Return list of linear devices in the netlist."""
        lin_devs = []
        for dev in self.devices:
            if not dev.is_nonlinear():
                lin_devs.append(dev)
        return lin_devs

    def get_nonlinear_devices(self):
        """Return list of nonlinear devices in the netlist."""
        nonlin_devs = []
        for dev in self.devices:
            if dev.is_nonlinear():
                nonlin_devs.append(dev)
        return nonlin_devs

    def get_node_idx(self, name):
        """
        Return a node index in the netlist from its name.

        Parameters
        ----------
        name : str
            Name of the node.

        Returns
        -------
        int
            Node index.

        """
        if name in self.node_name_to_idx:
            n = self.node_name_to_idx[name]
            return n
        else:
            logger.warning('Unknown node name: {}!'.format(name))
            return None

    def get_node_name(self, idx):
        """
        Return a node name from its index in the netlist.

        Parameters
        ----------
        idx : int
            Index of the node.

        Returns
        -------
        str
            Node name.

        """
        if idx < len(self.node_idx_to_name):
            n = self.node_idx_to_name[idx]
            return n
        else:
            logger.warning('Node index \'{}\' doesn\'t exist!'.format(idx))
            return None
    
    def get_mna_extra_rows_dict(self, analysis):
        """
        Return a mapping of the independent vsources to the device instances.

        Many devices require in their companion models additional rows to
        describe their behavior, associated with independent vsources. This
        method maps what are the MNA rows assigned to each device.

        Parameters
        ----------
        analysis : str
            Name of the analysis to get the dictionary from.

        Returns
        -------
        dict
            Maps an independent vsource index in the MNA formulation to the
            device instance responsible for it.

        """
        n = self.get_n()
        k = 0
        iidx = {}
        for dev in self.devices:
            if dev.get_num_vsources(analysis) > 0:
                iidx[dev] = n + k
                k = k + dev.get_num_vsources(analysis)
        return iidx

    def get_voltage_idx(self, name):
        """
        Return the row of the node voltage in the MNA formulation.

        Parameters
        ----------
        name : str
            Name of the node to get the MNA index.

        Returns
        -------
        int
            Node voltage index in the MNA formulation.

        """
        if name in self.node_name_to_idx:
            n = self.node_name_to_idx[name] - 1
            return n
        else:
            logger.warning('Unknown node name: {}!'.format(name))
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
        :class:`numpy.ndarray`
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
                if dev.get_num_vsources('dc') == 1:
                    print('{}:\t{:0.6f} A'.format(dev.name, x[n+k-1,0]))
                elif dev.get_num_vsources('dc') == 2:
                    print('{}:\t{:0.6f} A (I1)'.format(dev.name, x[n+k-1,0]))
                    print('{}:\t{:0.6f} A (I2)'.format(dev.name, x[n+k,0]))
                k = k + dev.get_num_vsources('dc')
            print()
        else:
            logger.warning('Unknown analysis name: {}!'.format(name))


    




