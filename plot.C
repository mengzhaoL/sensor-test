{ // simple plotting script for quick preview of the scanIV.py output
   TGraph *g = new TGraph("test.csv","%lf,%*lf,%lf");
   g->Draw("apl");
   g->SetMarkerStyle(8);
   g->GetXaxis()->SetTitle("Bias Voltage [V]");
   g->GetYaxis()->SetTitle("Measured Current [A]");
   gPad->SetLogy();
   g->GetYaxis()->SetRangeUser(1E-10,0.001);
}
