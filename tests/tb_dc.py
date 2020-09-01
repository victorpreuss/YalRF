import sys
sys.path.append("../yarf")

from yarf import Yarf

y = Yarf("Hello World!")

y.add_resistor('R1', 'n1', 'n2', 100)
y.add_resistor('R2', 'n2', 'gnd', 25)
y.add_resistor('R3', 'n3', 'n5', 50)
y.add_resistor('R4', 'n7', 'n8', 33)
y.add_resistor('R5', 'n8', 'n5', 66)
y.add_resistor('R6', 'n9', 'gnd', 50)

y.add_capacitor('C1', 'n4', 'n3', 1e-12)

y.add_inductor('L1', 'n4', 'n2', 1e-9)
y.add_inductor('L2', 'n5', 'gnd', 1e-9)

y.add_vdc('V1', 'n1', 'gnd', 10)
y.add_vdc('V2', 'n3', 'n2', 7)
y.add_idc('I1', 'n4', 'n1', 2)

y.add_vcvs('SRC1', 'n3', 'n7', 'n8', 'n5', 1.7)
y.add_ccvs('SRC2', 'n7', 'n9', 'gnd', 'n5', 12.4)

dc = y.add_dc_analysis('DC1')

y.run('DC1')

y.print_dc_voltages('DC1')
y.print_dc_currents('DC1')

