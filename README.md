# Sensor Test Scripts

## For users

### Installation

Install python3 and pyvisa
> pip3 install --user pyvisa

Download and install NI-VISA\
http://www.ni.com/zh-cn/support/downloads/drivers/download.ni-visa.html

Install a driver\
https://plugable.com/drivers/prolific/

Check out the sensor-test code.\
Connect the instruments to PC via USB.\
Check the VISA resource name in NI MAX.\
Update the VISA resource names in the .py files if needed.\
Call the testIO() function in python3 to check if your PC can talk to the instruments.\
For example, in a python3 session:
> import Kei2400CControl\
> smu=Kei2400CControl.keithley2400c()\
> smu.testIO()

### Perform a scan
An I-V scan
> python3 scanIV.py

### Visualize the results
For a quick preview of the test.csv file
> root -l plot.C

To superimpose several .csv files into a single plot, macros/plotall.py can be used.


## For developers 

- Fork the code with your personal github ID. See [details](https://help.github.com/articles/fork-a-repo/)

- Make your change, commit and push 

- Make a pull request. See [details](https://help.github.com/articles/using-pull-requests/)

## Some styles to follow 
- Minimize the number of main c++ files 
- Keep functions length less than one screen
- Seperate hard-coded cuts into script file 
- Document well the high-level bash file for work flow 
