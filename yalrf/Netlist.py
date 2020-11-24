from .Devices import *
from .Analyses import *
from .Utils import yalrf_logger as logger

import logging
logger.setLevel(logging.WARNING)

class Netlist():
    """
    Class containing the API for manipulating netlists.
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

    def add_vpulse(self, name, n1, n2, v1, v2, tstart, tstop=1e9, trise=1e-12, tfall=1e-12):
        """
        Add a pulse voltage source to the netlist. This source is used to apply
        voltage ramps at Transient simulations. It is 0 V for DC and AC analyses.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the voltage source.
        n2 : str
            Negative node of the voltage source.
        v1 : float
            Voltage level before tstart and after tstop.
        v2 : float
            Voltage level between tstart and tstop.
        tstart : float
            Time at which the voltage step is applied.
        tstop : float
            Time at which the voltage step ends.
        trise : float
            Rise time of the ramp at tstart.
        tfall : float
            Fall time of the ramp at tstop.

        Returns
        -------
        :class:`TransientVoltageSource`
            Reference to the created VoltageSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        vpulse = TransientVoltageSource(name, n1, n2, vtype='pulse', v1=v1, v2=v2, tstart=tstart, tstop=tstop, trise=trise, tfall=tfall)
        self.devices.append(vpulse)
        return vpulse

    def add_vsine(self, name, n1, n2, dc, ac, freq, phase=0):
        """
        Add a sinusoidal voltage source to the netlist. This source is used to apply
        sine waves at Transient simulations. It is 'dc' V for DC and 'ac' V for AC analyses.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Positive node of the voltage source.
        n2 : str
            Negative node of the voltage source.
        dc : float
            DC offset of the sine wave. Also used for DC bias point analysis.
        ac : float
            AC amplitude of the sine wave. Also used for AC analysis.
        freq : float
            Frequency of the sine wave at Transient analysis.
        phase : float
            Initial phase of the sine wave.

        Returns
        -------
        :class:`TransientVoltageSource`
            Reference to the created VoltageSource object.

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        vsine = TransientVoltageSource(name, n1, n2, vtype='sine', dc=dc, ac=ac, freq=freq, phase=phase)
        self.devices.append(vsine)
        return vsine

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

    def add_gyrator(self, name, n1, n2, n3, n4, G):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        n4 = self.add_node(n4)

        gyr = Gyrator(name, n1, n2, n3, n4, G)
        self.devices.append(gyr)
        return gyr

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
        the file Devices/BJT.py for a list of available parameters.

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

    def add_opamp(self, name, n1, n2, n3, G=100e3, Vmax=10e3):
        """
        Add an Operational Amplifier to the netlist.

        The instance of the added bjt is returned for the user. The opamp model
        is defined by a voltage gain and a maximum voltage for saturation.

        Parameters
        ----------
        name : str
            Name of the device.
        n1 : str
            Non-inverting input (+).
        n2 : str
            Inverting input (-).
        n3: str
            Output node.
        G: float
            Voltage gain of the opamp.
        Vmax: float
            Maximum voltage of the opamp. Above that, it saturates.

        Returns
        -------
        :class:`Opamp`
            Reference to the created Opamp object.

        Examples
        --------

        >>> opamp1 = y.add_opamp('Q1', 'n1', 'n2', 'n3', G=100, Vmax=20)

        """
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        
        opamp = Opamp(name, n1, n2, n3, G, Vmax)
        self.devices.append(opamp)
        return opamp

    # TODO:
    # def add_mosfet(self, name, n1, n2, n3, n4='gnd'):
    # def add_iprobe(self, name, n1, n2):
    # def add_transformer(self, name, n1, n2, n3, n4, T):
    # def add_subcircuit(self, name, YalRF object):
    # def remove_device(self, name):

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

    def get_num_nodes(self):
        """Return number of uniquely named nodes in the netlist."""
        return len(self.node_idx_to_name)

    def get_num_vsources(self):
        """
        Return number of independent voltage sources.

        The number of independent voltage sources needed are important for
        the MNA algorithm, since it determines the number of additional rows
        required in the matrix system. This number depends on the companion
        models used by each device.

        Returns
        -------
        int
            Number of independent voltage sources required.

        """
        m = 0
        for dev in self.devices:
            m = m + dev.get_num_vsources()
        return m

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
    
    def get_mna_extra_rows_dict(self):
        """
        Return a mapping of the independent vsources to the device instances.

        Many devices require in their companion models additional rows to
        describe their behavior, associated with independent vsources. This
        method maps what are the MNA rows assigned to each device.

        Returns
        -------
        dict
            Maps an independent vsource index in the MNA formulation to the
            device instance responsible for it.

        """
        n = self.get_num_nodes()
        k = 0
        iidx = {}
        for dev in self.devices:
            if dev.get_num_vsources() > 0:
                iidx[dev] = n + k
                k = k + dev.get_num_vsources()
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


