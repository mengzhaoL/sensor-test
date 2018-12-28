import ROOT
from ROOT import gPad,gStyle
from ROOT import TGraph,TCanvas,TPad,TMultiGraph

name = "EXX28995-WNo18_Single_Set-P6" # name of the sample set, for graph title and output file name
fname = [] # input csv file list
fname.append('LG1_SE2')
fname.append('LG1_SE2_2')
fname.append('LG1_SE3')
fname.append('LG1_SE3_2')
fname.append('LG1_SE5')
fname.append('LG1_SE5_2')
fname.append('LG1_SE5_3')
fname.append('LG1_SE5_4')
fname.append('LG1_SE5_NM')
fname.append('PIN1_SE5')
fname.append('PIN1_SE5_NM')


mg = TMultiGraph("mg","")

for iFile in fname:
   g = TGraph(iFile+".csv","%lf,%*lf,%lf");
   g.SetMarkerStyle(8)
   g.SetMarkerSize(0.8)
   mg.Add(g)
   print g,iFile

c = TCanvas("c", "canvas", 800, 600)
c.cd()
gPad.SetLogy()
mg.SetTitle(name+";Bias Voltage [V];Measured Current [A]") # set the graph title, x-axis title and y-axis title
mg.Draw("APL PMC PLC") # draw all the graphs using automatic coloring
mg.GetHistogram().GetYaxis().SetRangeUser(1E-10,0.001); # y-axis range
#mg.GetHistogram().GetXaxis().SetRangeUser(0,150); # x-axis range
c.BuildLegend(0.4,0.8-0.05*len(fname),0.6,0.8,"","PL") # generate the legend and define its position and style

c.Print(name+".pdf")
