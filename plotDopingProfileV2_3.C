//----------------------------------------------------
//               LGAD CV  plotDopingProfile  _IHEP
//
// Version 2.3   / Auto Get Depletion capacitance.   _mzli 20200801
// Version 2.2   / Process multiple files at once.   _mzli 20200729
// Version 2.1   / Add the smoothDoping() function & the TotalDoping() function.   _mzli 20200718
// Version 2     / Add  killPoints() function & smoothPoints() function.  _mzli 20200717
//----------------------------------------------------

#include <TH2.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>

TGraph *killPoints(TGraph *g, double startX, double endX, double startY, double endY); //Kill Range
TGraph *smoothPoints(TGraph *g);//  Five points   Fit
TGraph *smoothDoping(TGraph *g);//  Three-point mean
double TotalDoping(TGraph *g_dop,double D1,double D2);//  Calculate TotalDoping

TGraph *GetCVdata(TString dir,TString fname);//  1/C^2-V
TGraph *GetC_2(TGraph *g_cap);//  1/C^2-V
TGraph *GetDoping(TGraph *g_cap);//  doping profile
double GetCap(TGraph *g_Cap, TGraph *g_invcap2, double x1, double x2, double x3, double x4);


//TString dir = "/media/sf_Share2/LGAD/NDL/CV/B14/B14_3";  TString fname = "B14_3_mz_1_pad4";
//TString dir = "/media/sf_Share2/LGAD/NDL/CV/B14/B14_1";  TString fname = "B14_1_mz_1_pad1";
//TString dir2 = "/media/sf_Share2/LGAD/NDL/CV/B14/B14_1";  TString fname2 = "B14_1_mz_1_pad4";
//TString dir3 = "/media/sf_Share2/LGAD/NDL/CV/B14/B14_1";  TString fname3 = "B14_1_mz_1_pad3_pin";
//TString dir4 = "/media/sf_Share2/LGAD/NDL/CV/B14/B14_1";  TString fname4 = "B14_1_mz_1_pad4";
//TString dir = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_1";  TString fname = "B58_1_A2_pad1";
//TString dir2 = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_1";  TString fname2 = "B58_1_A2_pad2";
//TString dir3 = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_1";  TString fname3 = "B58_1_A2_pad3";
//TString dir4 = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_1";  TString fname4 = "B58_1_A2_pad4";
//TString dir = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_3";  TString fname = "B58_3_A4_pad4";
//TString dir2 = "/media/sf_Share2/LGAD/NDL/CV/B14/B14_1";  TString fname2 = "B14_1_mz_1_pad1";
//TString dir3 = "/media/sf_Share2/LGAD/NDL/CV/NDL9";  TString fname3 = "NDL9_mz_1_0721";

//TString dir = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_1";  TString fname = "B58_1_F5_pad1";
//TString dir2 = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_1";  TString fname2 = "B58_1_F5_pad2";
//TString dir3 = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_1";  TString fname3 = "B58_1_F5_pad4";
//TString dir4 = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_1";  TString fname4 = "B58_1_E3_pad4";
TString dir2 = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_1";  TString fname2 = "B58_1_A3_pad4";
TString dir = "/media/sf_Share2/LGAD/NDL/CV/For JSI irradiation/B58IVCV/B58_3_CV";  TString fname = "B58_3_G2_pad4_0V_good";
//TString dir3 = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_1";  TString fname3 = "B58_1_E3_pad4";
//TString dir = "/media/sf_Share2/LGAD/NDL/CV/B58/B58_3";  TString fname = "B58_3_A4_pad4";


double A = 0.13*0.13;//0.19*0.22; //  Sensor area [cm^2]
double esi = 11.7;
double e0 = 8.854E-14; // [F/cm]
double q = 1.602E-19; // [C]

void plotDopingProfileV2_3()
{

    //TCanvas *c0 = new TCanvas("c0","c0",2100,500);
    //c0->Divide (3,1);

    // -------------- Plot C-V
    TCanvas *c1 = new TCanvas("c1","c1",600,450);
    //c0->cd(1);
    TGraph *g_cap = GetCVdata(dir,fname); // CV data   (KIlled or SMOOTHed
    g_cap->SetMarkerStyle(24);
    g_cap->SetMarkerColor(4);
    g_cap->Draw("AP");
    gPad->SetLogy();
    TGraph *g_cap2 = GetCVdata(dir2,fname2); // CV data   (KIlled or SMOOTHed
    g_cap2->SetMarkerStyle(26);
    g_cap2->SetMarkerColor(2);
    g_cap2->SetLineWidth(0);
    g_cap2->Draw("CP");
    gPad->SetLogy();
    /*TGraph *g_cap3 = GetCVdata(dir3,fname3); // CV data   (KIlled or SMOOTHed
    g_cap3->SetMarkerStyle(26);
    g_cap3->SetMarkerColor(8);
    g_cap3->SetLineWidth(0);
    g_cap3->Draw("CP");
    gPad->SetLogy();
    TGraph *g_cap4 = GetCVdata(dir4,fname4); // CV data   (KIlled or SMOOTHed
    g_cap4->SetMarkerStyle(27);
    g_cap4->SetMarkerColor(46);
    g_cap4->SetLineWidth(0);
    g_cap4->Draw("CP");
    gPad->SetLogy(); //*/

    TLegend *legendC = new TLegend(0.60,0.80,0.81,0.90);
    legendC->AddEntry(g_cap,"50 microns, B58-3","p");
    legendC->AddEntry(g_cap2,"33 microns, B58-1","p");
    //legendC->AddEntry(g_cap3,fname3,"p");
    //legendC->AddEntry(g_cap4,fname4,"p");
    legendC->Draw();

    // -------------- Plot 1/C^2-V
    TCanvas *c2 = new TCanvas("c2","c2",600,450);
    //c0->cd(2);
    TGraph *g_invcap2 = GetC_2(g_cap);//  ------------------  1/C^2-V
    g_invcap2->SetMarkerStyle(24);
    g_invcap2->SetMarkerColor(4);
    g_invcap2->Draw("AP");
    //cout<<GetCap(g_cap,g_invcap2, 40, 70, 85, 100)<<" <<-- sensor1 dCap "<<endl;
    TGraph *g_invcap22 = GetC_2(g_cap2);//  ------------------ 1/C^2-V
    g_invcap22->SetMarkerStyle(26);
    g_invcap22->SetMarkerColor(2);
    g_invcap22->SetLineWidth(0);
    g_invcap22->Draw("CP");
    //cout<<GetCap(g_cap2,g_invcap22, 30, 40, 50, 65)<<" <<-- sensor2 dCap "<<endl;
    /*TGraph *g_invcap23 = GetC_2(g_cap3);//  -------------------  1/C^2-V
    g_invcap23->SetMarkerStyle(26);
    g_invcap23->SetMarkerColor(8);
    g_invcap23->SetLineWidth(0);
    g_invcap23->Draw("CP");
    //cout<<GetCap(g_cap3,g_invcap23, 35, 45, 50, 70)<<" <<-- sensor3 dCap "<<endl;
    //TGraph *g_invcap24 = GetC_2(g_cap4);//  -------------------  1/C^2-V
    //g_invcap24->SetMarkerStyle(27);
    //g_invcap24->SetMarkerColor(+46);
    //g_invcap24->SetLineWidth(0);
    //g_invcap24->Draw("CP");
    //cout<<GetCap(g_cap4,g_invcap24, 35, 45, 50, 70)<<" <<-- sensor4 dCap "<<endl;//*/

    TLegend *legendC2 = new TLegend(0.20,0.80,0.41,0.90);
    legendC2->AddEntry(g_invcap2,"50 microns, B58-3","p");
    legendC2->AddEntry(g_invcap22,"33 microns, B58-1","p");
    //legendC2->AddEntry(g_invcap23,fname3,"p");
    //legendC2->AddEntry(g_invcap24,fname4,"p");
    legendC2->Draw();

    // -------------- Plot doping profile
    TCanvas *c3 = new TCanvas("c3","c3",600,450);
    //c0->cd(3);
    TGraph *g_dop = GetDoping(g_cap);//  doping profile
    //g_dop = smoothDoping(g_dop);   //smooth Doping
    g_dop->SetMarkerStyle(24);
    g_dop->SetMarkerColor(4);
    g_dop->Draw("AP");
    //gPad->SetLogy();
    TGraph *g_dop2 = GetDoping(g_cap2);//  doping profile
    g_dop2->SetMarkerStyle(26);
    g_dop2->SetMarkerColor(2);
    g_dop2->SetLineWidth(0);
    g_dop2->Draw("CP");
    //gPad->SetLogy();
    /*TGraph *g_dop3 = GetDoping(g_cap3);//  doping profile
    g_dop3->SetMarkerStyle(26);
    g_dop3->SetMarkerColor(8);
    g_dop3->SetLineWidth(0);
    g_dop3->Draw("CP");
    //TGraph *g_dop4 = GetDoping(g_cap4);//  doping profile
    //g_dop4->SetMarkerStyle(27);
    //g_dop4->SetMarkerColor(46);
    //g_dop4->SetLineWidth(0);
    //g_dop4->Draw("CP"); //*/

    TLegend *legendD = new TLegend(0.60,0.80,0.81,0.90);
    legendD->AddEntry(g_dop,"50 microns, B58-3","p");
    legendD->AddEntry(g_dop2,"33 microns, B58-1","p");
    //legendD->AddEntry(g_dop3,fname3,"p");
    //legendD->AddEntry(g_dop4,fname4,"p");
    legendD->Draw();


    c1->SaveAs("C_V.pdf","recreate");
    c2->SaveAs("C2_V.pdf");
    c3->SaveAs("Doping.pdf");




/*
     //  Kill  Points
    TGraph *gKPed =killPoints(g_capTemp,-1, 100,5000,110E+48);  //Kill Range
    //gKPed =killPoints(gKPed,20, 25,0,90);    //Kill again

    gKPed->SetTitle("killPoints");
    c0->cd(2);
    gKPed->Draw("AC*");

     //  smoothPoints
    TGraph *gSPed = smoothPoints(gKPed);
    //gSPed = smoothPoints(gSPed);     //again
    //gSPed = smoothPoints(gSPed);
    gSPed->SetTitle("smoothPoints");
    c0->cd(3);
    gSPed->Draw("AC*");
*/

}



//-----------------------------------------------My Function---------------------------------
TGraph *GetCVdata(TString dir,TString fname){
    TGraph *g_capTemp = new TGraph( Form("%s/%s.csv",dir.Data(),fname.Data()), "%lf,%*lf,%lf"); //read Raw data
    int N = g_capTemp->GetN();
    double *volTemp = g_capTemp->GetX(); // [V]
    double *capTemp = g_capTemp->GetY(); // [pF]
    for(int i=0; i<N; ++i){
        volTemp[i] = -volTemp[i]; // convert to positive
        capTemp[i] = capTemp[i]*1E-12;   //   [pF]->[F]
        g_capTemp->SetPoint(i, volTemp[i], capTemp[i]);
    }

    //g_capTemp = smoothPoints(g_capTemp);

    g_capTemp->SetTitle("C_V;Bias Voltage [V];Capacitance [F]" );

    return g_capTemp;
}

TGraph *GetC_2(TGraph *g_cap){

    int nvol = g_cap->GetN();
    double *vol = g_cap->GetX();
    double *cap = g_cap->GetY();

    double *invcap2 = new double[nvol]; // 1/C^2
   // double *invcap2_dv = new double[nvol-1]; // d(1/C^2)/dV

    for(int i=0; i<nvol; ++i){
        invcap2[i] = 1.0/cap[i]/cap[i]; // [[1/F^2]
    }

    TGraph *g_invcap2 = new TGraph(nvol, vol, invcap2);
    //TGraph *g_invcap2_dv = new TGraph(nvol-1, vol, invcap2_dv);
    g_invcap2->SetTitle( "1/C^2;Bias Voltage [V];Capacitance^{-2} [F^{-2}]" );
    //g_invcap2->Draw("AP*");
    //g_invcap2->SetLogx();
    //gPad->SetLogy();

    return g_invcap2;
}

double GetCap(TGraph *g_Cap, TGraph *g_invcap2, double x1, double x2, double x3, double x4){

    double f1p0,f1p1,f2p0,f2p1;
    double dVol,dCap;
    TF1 *f1=new TF1("f1","[0]+[1]*x",x1,x2);
    f1->SetParameters(-3E+22,1E+21);
    g_invcap2->Fit(f1,"R");
    double *par1 = f1->GetParameters();
    f1p0 = par1[0]; f1p1 = par1[1];\
    cout<<"f1  p  "<<par1[0]<<" "<<par1[1]<<endl;

    //TF1 *f2=new TF1("f2","[0]+[1]*x",85,95);
    TF1 *f2=new TF1("f2","pol1",x3,x4);
    //f1->SetParameters(2.2E+22,2.8E+19);
    g_invcap2->Fit(f2,"R+");
    double *par2 = f2->GetParameters();
    f2p0 = par2[0]; f2p1 = par2[1];
    cout<<"f2  p  "<<par2[0]<<" "<<par2[1]<<endl;

    dVol = (f2p0-f1p0)/(f1p1-f2p1);

    //TF1 *f3 = new TF1("f3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", dVol-5,dVol+5);
    TF1 *f3 = new TF1("f3", "pol3", dVol-5,dVol+5);
    g_Cap->Fit(f3,"R+");
    double *par3 = f3->GetParameters();
    //if(f->GetChisquare()>1e-4) fitStatus=1;
    dCap = par3[0]+par3[1]*dVol+par3[2]*dVol*dVol+par3[3]*dVol*dVol*dVol;
     cout<<" DV DC "<<dVol<<" "<<dCap<<endl;
    return dCap;
}


TGraph *GetDoping(TGraph *g_cap){

    int nvol = g_cap->GetN();
    double *vol = g_cap->GetX();
    double *cap = g_cap->GetY();

    double *invcap2 = new double[nvol]; // 1/C^2
    double *invcap2_dv = new double[nvol-1]; // d(1/C^2)/dV
    double *doppf = new double[nvol-1];
    double *depth = new double[nvol-1];

    for(int i=0; i<nvol; ++i){
        invcap2[i] = 1.0/cap[i]/cap[i]; // [[1/F^2]
    }

    // Calculate derivatives and doping profile
    for(int i=0; i<nvol-1; ++i){
        invcap2_dv[i] = fabs(invcap2[i+1]-invcap2[i])/(vol[i+1]-vol[i]); // forward difference
        doppf[i] = 2.0/q/esi/e0/A/A/invcap2_dv[i]; // [cm^{-3}]
        //if (doppf[i]> 50E+15) doppf[i] = doppf[i-1];  //kill bad points
        depth[i] = A*esi*e0*(1.0/cap[i])*1E+4; // [cm] to [um]
        //if(depth[i]<0) depth[i]=depth[i-1];
    }

    TGraph *g_dop = new TGraph(nvol-1, depth, doppf);
    g_dop = smoothDoping(g_dop);   //smooth Doping

    //double totalDoping = TotalDoping(g_dop, 1,  3); // Fit and Set Integration range [um]
    //cout<<" ^_^  -  ^_^  -  ^_^   TotalDoping  "<<totalDoping*1E-4<<" cm^{-3} "<<endl; //  um->cm

    g_dop->SetTitle("Doping profile;Depth [um];Doping profile [cm^{-3}]");
    g_dop->GetXaxis()->SetRangeUser(0,3);

    return g_dop;

}

TGraph *killPoints(TGraph *g,double startX, double endX,double startY, double endY){

    int nvol = g->GetN();
    double *vol = g->GetX(); // [V]
    double *cap = g->GetY(); // [pF]

    TGraph *gKilled = new TGraph();

    int newCount=0;
    for(int i =0; i<nvol; i++){
        if(!(vol[i]>startX && vol[i]<endX  && cap[i]>startY && cap[i]<endY)){
            gKilled->SetPoint(newCount,vol[i],cap[i]);
            newCount++;
        }
    }

    //TCanvas *cK = new TCanvas("cK","cK",1000,500);

    //gKilled->Draw("AC*");

    return gKilled;
    //delete gKilled;
}

TGraph *smoothPoints(TGraph *g){

    int nvol = g->GetN();
    double *vol = g->GetX();
    double *cap = g->GetY();

    TGraph *gSmoothed = new TGraph();

    double  C[5],V[5];

    for(int i=0; i<nvol; ++i){

        if(i>1 && i<nvol-2){
            C[0] = cap[i-2]; V[0] = vol[i-2];
            C[1] = cap[i-1]; V[1] = vol[i-1];
            C[2] = cap[i];   V[2] = vol[i];
            C[3] = cap[i+1]; V[3] = vol[i+1];
            C[4] = cap[i+2]; V[4] = vol[i+2];
            //cout<<"-----------------cap "<<C[0]<<"   "<<V[0]<<endl;
            //cout<<"-----------------cap "<<C[1]<<"   "<<V[1]<<endl;
            //cout<<"-----------------cap "<<C[2]<<"   "<<V[2]<<endl;
            //cout<<"-----------------cap "<<C[3]<<"   "<<V[3]<<endl;
            //cout<<"-----------------cap "<<C[4]<<"   "<<V[4]<<endl;

            //TCanvas *cS = new TCanvas("cS","cS",800,600);
            TGraph *g=new TGraph(5,V,C);
            //TF1 *f=new TF1("f","pol2",V[0]-0.01,V[4]+0.01);
            TF1 *f = new TF1("f", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", V[0]-0.01,V[4]+0.01);
            //g->Draw("AC*");//"AL*");
            g->Fit(f,"R+");
            double *par = f->GetParameters();
            //if(f->GetChisquare()>1e-4) fitStatus=1;
            double cap_i = par[0]+par[1]*V[2]+par[2]*V[2]*V[2]+par[3]*V[2]*V[2]*V[2];
            gSmoothed->SetPoint(i,vol[i],cap_i);
        }
        else gSmoothed->SetPoint(i,vol[i],cap[i]);

    }
    //TCanvas *cSS = new TCanvas("cSS","cSS",1000,500);
    //gSmoothed->Draw("AC*");
     return gSmoothed;
}

TGraph *smoothDoping(TGraph *g){

    int n = g->GetN();
    double *Deep = g->GetX();
    double *Doping = g->GetY();

    TGraph *gSmoothed = new TGraph();

    for(int i=0; i<n; ++i){
        if (i >0 && i < n-2) {
            double newdoppf_i = (Doping[i-1]+Doping[i]+Doping[i+1])/3; //
            gSmoothed->SetPoint(i,Deep[i],newdoppf_i);
        }
        else  gSmoothed->SetPoint(i,Deep[i],Doping[i]);;
    }
    return gSmoothed;
}


double TotalDoping(TGraph *g_dop, double D1, double D2){

    TF1 *f4=new TF1("f4","gaus",D1,D2);
    g_dop->Fit(f4,"R+");
    double Integral = f4->Integral(D1,D2, 1.e-12);

    return Integral;
}

