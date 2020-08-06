from .Devices import *
from .Analyses import *

class Yarf():
    
    def __init__(self, name):
        self.name = name # circuit name given by the user

        # netlist related parameters
        self.devices = list()   # list containing all the devices in the circuit
        self.nodes = dict()     # dictionary to associate node names to index
        self.node_name = list() # list of all the node names
        self.nodes['gnd'] = 0   # initialize ground node
        self.node_name.append('gnd')

        # analyses related parameters
        self.analyses = dict() # dictionary of analysis to run in the circuit

    def add_node(self, n):
        if n not in self.nodes:
            self.nodes[n] = len(self.nodes)
            self.node_name.append(n)

        return self.nodes[n]
        
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

    def add_vdc(self, name, n1, n2, dc):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        vdc = VoltageSource(name, n1, n2, vtype='dc', dc=dc)
        self.devices.append(vdc)
        return vdc

    def add_idc(self, name, n1, n2, dc):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        idc = CurrentSource(name, n1, n2, itype='dc', dc=dc)
        self.devices.append(idc)
        return idc

    def add_vsource(self, name, n1, n2, vtype=None, dc=None, mag=None, phase=None):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        vsource = VoltageSource(name, n1, n2, vtype, dc, mag, phase)
        self.devices.append(vsource)
        return vsource

    def add_isource(self, name, n1, n2, itype=None, dc=None, mag=None, phase=None):
        n1 = self.add_node(n1)
        n2 = self.add_node(n2)
        
        isource = CurrentSource(name, n1, n2, itype, dc, mag, phase)
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

    # TODO: implement this method using a dictionary for the devices
    def get_device(self, name):
        print('WARNING: method not implemented!')
        return None

    def get_analysis(self, name):
        if name in self.analyses:
            return self.analyses[name]
        else:
            print('WARNING: analysis name \'{}\' already taken!'.format(name))
            return None

    def add_dc_analysis(self, name):
        if name not in self.analyses:
            self.analyses[name] = DC(name)
            return self.analyses[name]
        else:
            print('WARNING: analysis name \'{}\' already taken!'.format(name))
            return None

    def print_dc_voltages(self, name):
        if name in self.analyses:
            self.analyses[name].print_dc_voltages(self)
        else:
            print('WARNING: unknown analysis name: {}'.format(name))
    
    def print_dc_currents(self, name):
        if name in self.analyses:
            self.analyses[name].print_dc_currents(self)
        else:
            print('WARNING: unknown analysis name: {}'.format(name))

    def run(self, name):
        if name in self.analyses:
            return self.analyses[name].run(self)
        else:
            print('WARNING: unknown analysis name: {}'.format(name))

    def run_all(self):
        for a in self.analyses.items():
            a.run(self)

    def is_linear(self):
        for d in self.devices:
            if d.is_linear() == False:
                return False
        return True




