import KeyE4980AControl
import Kei2400CControl as kei2400
import csv
import numpy as np
import platform
import sys
import warnings
import time

def csv_writer(data,path):
    with open(path,"w") as csv_file:
        writer=csv.writer(csv_file,lineterminator='\n')
        writer.writerow(['Bias Voltage[V]','Measured Voltage[V]','Measured Capacitance[pF]','Capacitance^{-2}[pF]^{-2}'])
        for val in data:
            writer.writerows([val])

# prevent from running with python2
if platform.python_version().startswith('2'):
   print('You are running with',platform.python_version(),'.')
   print('Please run with python3.')
   print('Exit.')
   sys.exit()

### LCR meter settings
lcr=KeyE4980AControl.keysighte4980a("USB0::0x2A8D::0x2F01::MY46516486::INSTR")
lcr.set_frequency("1MHz")
lcr.set_voltage_level(0.1) # in V
timedelay=1.0 # in second

### Source meter settings
biasSupply=kei2400.keithley2400c("ASRL1::INSTR")
biasSupply.set_current_protection(100E-6) # current protection in A
biasSupply.set_voltage_protection(200) # voltage protection in V
positiveHV=True # sign of the voltage
HVrange=80.0*1e3  # voltage scan range in mV in absolute value


vols=[]
mvols=[]
pcap=[]
invpcap2=[]

if positiveHV:
    sign=1
else:
    sign=-1

iStart=int(0*1e3)
iEnd=int(sign*HVrange+sign*1)
iStep=int(sign*1.0*1e3)
lcr.set_trigger_remote()
biasSupply.output_on()
for iBias in range(iStart,iEnd,iStep):
    biasvol=iBias/1000 # mV to V
    if biasvol>200:
        warnings.warn("Warning! Exceed the maximum voltage of the bias adapter!")
        break
    vols.append(biasvol)
    mvols.append(biasSupply.set_voltage(biasvol))
    time.sleep(timedelay)
    ipcap=lcr.get_capacitance()/1e-12
    pcap.append(ipcap)
    invpcap2.append(1.0/(ipcap*ipcap))
    if biasSupply.hit_compliance():
        break

print("Bias Vols: "+str(vols))
print("Measured vols: "+str(mvols))
print("Capacitance: "+str(pcap))

data=[vols,mvols,pcap,invpcap2]
dataarray=np.array(data)

filename="test.csv"
csv_writer(dataarray.T,filename)

lcr.set_trigger_internal()
biasSupply.set_voltage(0*1e3)
biasSupply.output_off()
