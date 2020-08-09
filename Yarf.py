from .Devices import *
from .Analyses import *

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('[%(levelname)s]: %(name)s: %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

logger.addHandler(stream_handler)

class Yarf():
    
    def __init__(self, name):
        self.name = name # circuit name given by the user

        # netlist related parameters
        self.devices = []           # list of all the devices
        self.node_idx = {'gnd': 0}  # node names to index
        self.node_name = ['gnd']    # index to node name

        # dictionary of analysis to run in the circuit
        self.analyses = {} 

    def get_devices(self):
        return self.devices

    def get_device(self, name):
        for dev in self.devices:
            if dev.name == name:
                return dev
        return None

    def get_v_idx(self, name):
        return (self.node_idx[name] - 1)

    def get_analysis(self, name):
        if name in self.analyses:
            return self.analyses[name]
        else:
            logger.warning('Unknown analysis name: {}!'.format(name))
            return None

    def add_ac_analysis(self, name, start, stop, numpts=10, stepsize=None, sweeptype='linear'):
        if name not in self.analyses:
            self.analyses[name] = AC(name, start, stop, numpts, stepsize, sweeptype)
            return self.analyses[name]
        else:
            logger.warning('Analysis name \'{}\' already taken!'.format(name))
            return None

    def add_dc_analysis(self, name):
        if name not in self.analyses:
            self.analyses[name] = DC(name)
            return self.analyses[name]
        else:
            logger.warning('Analysis name \'{}\' already taken!'.format(name))
            return None

    def run(self, name, x0=None):
        if name in self.analyses:
            return self.analyses[name].run(self, x0)
        else:
            logger.warning('Unknown analysis name: {}!'.format(name))
            return None

    def add_node(self, n):
        if n not in self.node_idx:
            self.node_idx[n] = len(self.node_idx)
            self.node_name.append(n)
        return self.node_idx[n]
        
    def add_resistor(self, name, n1, n2, value):
        n1 = self.add_node(n1) # add node string and get the
        n2 = self.add_node(n2) # equivalent index for the MNA
        
        res = Resistor(name, n1, n2, value)
        self.devices.append(res)
        return res
        
    def add_capacitor(self, name, n1, n2, value):  
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        cap = Capacitor(name, n1, n2, value)
        self.devices.append(cap)
        return cap
        
    def add_inductor(self, name, n1, n2, value):  
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        ind = Inductor(name, n1, n2, value)
        self.devices.append(ind)
        return ind

    def add_vdc(self, name, n1, n2, dc=0):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        vdc = VoltageSource(name, n1, n2, vtype='dc', dc=dc)
        self.devices.append(vdc)
        return vdc

    def add_idc(self, name, n1, n2, dc=0):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        idc = CurrentSource(name, n1, n2, itype='dc', dc=dc)
        self.devices.append(idc)
        return idc

    def add_vsource(self, name, n1, n2, vtype=None, dc=0, ac=0, phase=0):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        vsource = VoltageSource(name, n1, n2, vtype, dc, ac, phase)
        self.devices.append(vsource)
        return vsource

    def add_isource(self, name, n1, n2, itype=None, dc=0, ac=0, phase=0):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        isource = CurrentSource(name, n1, n2, itype, dc, ac, phase)
        self.devices.append(isource)
        return isource

    def add_vcvs(self, name, n1, n2, n3, n4, G=1.0, tau=None):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        n4 = self.add_node(n4)

        vcvs = VoltageControlledVoltageSource(name, n1, n2, n3, n4, G, tau)
        self.devices.append(vcvs)
        return vcvs

    def add_vccs(self, name, n1, n2, n3, n4, G=1.0, tau=None):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        n4 = self.add_node(n4)

        vccs = VoltageControlledCurrentSource(name, n1, n2, n3, n4, G, tau)
        self.devices.append(vccs)
        return vccs
        
    def add_ccvs(self, name, n1, n2, n3, n4, G=1.0, tau=None):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        n4 = self.add_node(n4)

        ccvs = CurrentControlledVoltageSource(name, n1, n2, n3, n4, G, tau)
        self.devices.append(ccvs)
        return ccvs
        
    def add_cccs(self, name, n1, n2, n3, n4, G=1.0, tau=None):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        n4 = self.add_node(n4)

        cccs = CurrentControlledCurrentSource(name, n1, n2, n3, n4, G, tau)
        self.devices.append(cccs)
        return cccs

    def add_diode(self, name, n1, n2):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        diode = Diode(name, n1, n2)
        self.devices.append(diode)
        return diode

    def add_bjt(self, name, n1, n2, n3):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        n3 = self.add_node(n3)
        
        bjt = BJT(name, n1, n2, n3)
        self.devices.append(bjt)
        return bjt

    # returns the number of uniquely named nodes in the circuit
    def get_n(self):
        return len(self.node_name)

    # returns the number of independent voltages
    def get_m(self, analysis):
        m = 0
        for dev in self.devices:
            m = m + dev.get_num_vsources(analysis)
        return m

    # given a node name, returns its index in the netlist
    def get_node_idx(self, name):
        if name in self.node_idx:
            return self.node_idx[name]
        else:
            return None

    # given a node index in the netlist, returns its name
    def get_node_name(self, idx):
        if idx < len(self.node_name):
            return self.node_name[idx]
        else:
            return None

    # return a list of the linear devices in the netlist
    def get_linear_devices(self):
        lin_devs = []
        for dev in self.devices:
            if dev.is_linear():
                lin_devs.append(dev)
        return lin_devs

    # return a list of the nonlinear devices in the netlist
    def get_nonlinear_devices(self):
        nonlin_devs = []
        for dev in self.devices:
            if not dev.is_linear():
                nonlin_devs.append(dev)
        return nonlin_devs

    # return True if netlist is purely linear
    def is_linear(self):
        for dev in self.devices:
            if dev.is_linear() == False:
                return False
        return True
    
    # return a dict which maps indep vsource idx in the MNA to a device
    def get_mna_extra_rows_dict(self, analysis):
        n = self.get_n()
        k = 0
        iidx = {}
        for dev in self.devices:
            if dev.get_num_vsources(analysis) > 0:
                iidx[dev] = n + k
                k = k + dev.get_num_vsources(analysis)
        return iidx

    def print_dc_voltages(self, name):
        if name in self.analyses:
            x = self.analyses[name].get_dc_solution()
            for i in range(1, self.get_n()):
                print('{}:\t{:0.4f} V'.format(self.node_name[i], x[i-1,0]))
            print()
        else:
            logger.warning('Unknown analysis name: {}!'.format(name))
    
    def print_dc_currents(self, name):
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



