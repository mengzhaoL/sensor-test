# Sensor Test Scripts

## For users

### Installation

Install python3, pyvisa and matplotlib
> pip3 install --user pyvisa \
> pip3 install --user matplotlib

Download and install NI-VISA\
http://www.ni.com/zh-cn/support/downloads/drivers/download.ni-visa.html

Install a driver\
https://plugable.com/drivers/prolific/

Check out the sensor-test code.\
Connect the instruments to PC via USB.\
Check the VISA resource name in NI MAX.\
Update the VISA resource names in the .py files if needed.\
Check if your PC can talk to the instruments.\
For example,
> python3 Kei2400CControl.py

should print out the identification information like KEITHLEY INSTRUMENTS...

### Perform a scan
An I-V scan with a single source meter:
> python3 scanIV.py

An I-V scan with two source meters (Keithley 2410 as power supply and Keithley 2400 as current meter):\
(Note: Mac OS is not supported)
> python3 scanIV2.py

A C-V scan:
> python3 scanCV.py

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
