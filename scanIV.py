import Kei2400CControl as kei2400
import visa
import time
import pylab
import csv
import numpy as np
import platform
import sys

def csv_writer(data,path):
    with open(path,"w") as csv_file:
        writer=csv.writer(csv_file,lineterminator='\n')
        writer.writerow(['Bias Voltage[V]','Measure Voltage [V]','Meaure Current [A]'])
        for val in data:
            writer.writerows([val])

# prevent from running with python2
if platform.python_version().startswith('2'):
   print('You are running with',platform.python_version(),'.')
   print('Please run with python3.')
   print('Exit.')
   sys.exit()

biasSupply=kei2400.keithley2400c()

vols=[]
mvols=[]
current=[]

iStart=int(0*1e3)
iEnd=int(5.0*1e3)
iStep=int(0.5*1e3)
for iBias in range(iStart,iEnd,iStep):
    biasSupply.output_on()
    biasvol=iBias/1000
	#if biasvol>2:
    #    break
    vols.append(biasvol)
    mvols.append(biasSupply.set_voltage(biasvol))
    current.append(biasSupply.display_current())

print("Bias Vols: "+str(vols))
print("Measure vols: "+str(mvols))
print("Current: "+str(current))

data=[vols,mvols,current]
dataarray=np.array(data)

filename="test.csv"
csv_writer(dataarray.T,filename)

biasSupply.set_voltage(0*1e3)
