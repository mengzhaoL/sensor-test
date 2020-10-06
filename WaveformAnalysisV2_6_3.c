//----------------------------------------------------
//               Waveform Analysis _mzli_IHEP
//
// Version 2.6.3   / Optimized the writing of histogram.
// Version 2.6.2   / Optimized the Charge_distribution() function.
// Version 2.6.1   / Save the canvas.
// Version 2.6   / Add  Charge_distribution() function.
// Version 2.5   / Add  Noise_distribution() function.
// Version 2.4.4 / Add noise for WF2 waveform in get_CFDt() function.
// Version 2.4.3 / The T_CFD1() function is optimized.
// Version 2.4.2 / Add the .csv file reading, time resolution get by CFDt1-CFDt2.
//----------------------------------------------------
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include<algorithm>
#include <iostream>

 //--
#define STEP 402                    //Step of waveform
#define DELTA_TIME 0.025             //Interval between two point [ns]
#define THRESHOLD 15               //Threshold [mV]
#define N 10000                      //Entries of waveform
#define BASELINE_START 0           //Range of baseline [ns]
#define BASELINE_STOP 4
#define Q_THRESHOLD 2*200 //000             //Threshold(By Q) [fC]
#define GATE_WIDHT 1000              //point number   Gate  = number*DELTA_TIME   */

/*//Just for WF2
#define STEP 50000                    //Step of waveform
#define DELTA_TIME 0.0001             //Interval between two point [ns]
#define THRESHOLD 1                //Threshold [mV]
#define N 1000//1000                     //Entries of waveform
#define BASELINE_START 0.001           //Range of baseline [ns]
#define BASELINE_STOP 0.003
#define Q_THRESHOLD 2000             //Threshold(By Q) [fC]
#define GATE_WIDHT 1000              //point number   Gate  = number*DELTA_TIME    */

//string title="/media/sf_Share2/LGAD/IME/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_180V/C2BiasV";             //file
//string title2="/media/sf_Share2/LGAD/IME/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_180V/C3BiasV";             //file
//string title="/media/sf_Share2/LGAD/IME/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_120V/C2BiasV";             //file
//string title2="/media/sf_Share2/LGAD/IME/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_120V/C3BiasV";             //file
//string title="/media/sf_Share2/LGAD/IME/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_140V/C2BiasV";             //file
//string title2="/media/sf_Share2/LGAD/IME/B58_1_A3_W8-IV-E4-L1-15_70_beta/200V_140V/C3BiasV";             //file
string title="/media/sf_Share2/LGAD/IME/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_200V/C2BiasV";             //file
string title2="/media/sf_Share2/LGAD/IME/B58_1_A3_W7-II-F5-L1-15_70_beta/220V_200V/C3BiasV";             //file


double *get_peak_parameters(double *waveform_point,double *peak_fit,int threshold,int step);
double get_baseline(double *waveform_point, double win_start, double win_stop);
//double get_tot(TString file,double threshold,int step);
double get_charge(TString file,double baseline,int step,int gate_width);    //Q within the gate width
//int get_rise_time_point_num(TString file,double min_v,double max_v,int step);
void get_waveform_data(TString file,double *waveform,int step);
double get_CFDt1(double *waveform,double v_baseline,double CFD,double PeakV,double PeakP);//CFD time
double get_CFDt2(double *waveform,double v_baseline,double CFD,double PeakV,double PeakP);//CFD time

//void TTS1();
void T_CFD1();
void T_CFD2();
//void TTS3();
//void ADC();
//void QDC();
void Baseline_distribution();
void Charge_distribution(double gate_width);//Q within the gate width
void Peak_distribution();
void Noise_distribution(double win_start,double win_stop);

void WaveformAnalysisV2_6_3(){
        cout<<" * * * * * * * * * * * * * * * * * * * Check0 * OK * * * * * * * * * * * * *"<<endl;
	//TTS1();                                  //TTS by peak time
        //T_CFD1();                                   //Time resolution by CFD  (t1) //For WF2 output files
        //T_CFD2();                                    //Time resolution by CFD  (t1-t2)
	//TTS3();                                  //TTS by T-A corrention
	//ADC();                                   //ADC Spectrum
	//Baseline_distribution();                 //Baseline
        Charge_distribution(4);                   // Charge_distribution      Q within the gate width [ns]
        //Peak_distribution();                      //peak spectrum
        //Noise_distribution(0.1,2.1);                //Noise spectrum.  The range [ns]

        //double waveform_point1[STEP];
        //TString file1="BV170_130V_10Wfm_Ch1";
        //get_waveform_data(file1,waveform_point1,STEP);
}

void T_CFD1(){

        TString file1;
        double CFD = 0.1;                    //CFD Fraction [0-1]
        int FillCounter =0;
        int j;
        double *peak1,PeakV1,PeakP1;
        double CFDt1,CFDt,peak_fit1[2],v_baseline1,q;
        static double waveform_point1[STEP];
        TH1F *h1=new TH1F("Cap10_0.1","CFD-0.1",500,-1,1);//0-1ns
        TH1F *h2=new TH1F("Cap10_0.2","CFD-0.2",500,-1,1);//0-1ns
        TH1F *h3=new TH1F("Cap10_0.3","CFD-0.3",500,-1,1);//0-1ns
        TH1F *h4=new TH1F("Cap10_0.4","CFD-0.4",500,-1,1);//0-1ns
        TH1F *h5=new TH1F("Cap10_0.5","CFD-0.5",500,-1,1);//0-1ns
        TH1F *h6=new TH1F("Cap10_0.6","CFD-0.6",500,-1,1);//0-1ns
        TH1F *h7=new TH1F("Cap10_0.7","CFD-0.7",500,-1,1);//0-1ns
        TH1F *h8=new TH1F("Cap10_0.8","CFD-0.8",500,-1,1);//0-1ns
        TH1F *h9=new TH1F("Cap10_0.9","CFD-0.9",500,-1,1);//0-1ns
        for(j=0;j<N;j++){
             file1=title+to_string(j+1);    //
             cout<<" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  B e g i n  "<<j+1<<endl;
             get_waveform_data(file1,waveform_point1,STEP);
             //cout<<" - - - - - - - - - - - - - - - - - - - - - - - - - check-1"<<endl;
             peak1 = get_peak_parameters(waveform_point1,peak_fit1,THRESHOLD,STEP);
             PeakV1 = peak1[0];  PeakP1 = peak1[1];
             cout<<"--PeakV1 "<<PeakV1<<"mV"<<endl;
             //cout<<"--PeakFit  :"<<peak_fit1[0]<<" mV"<<endl;
             v_baseline1 = get_baseline(waveform_point1,BASELINE_START,BASELINE_STOP);

             if(PeakV1>1*THRESHOLD && PeakV1<150*THRESHOLD){
                   FillCounter++;
                   cout<<"---------------------------------------------------------------------- Fill counter "<<FillCounter<<endl;
                  for(CFD = 0.1;CFD<0.91;CFD+=0.1){
                     //CFD = 0.8;
                     CFDt1 = get_CFDt1(waveform_point1,v_baseline1,CFD,PeakV1,PeakP1);
                     CFDt = CFDt1;

                     int intCFD = CFD*10+0.1;
                     cout<<" - - - intCFD "<<intCFD<<"  CFD "<<CFD<<endl;

                     switch (intCFD)
                     {
                         case 1 :
                             h1->Fill(CFDt);;
                             break;
                         case 2 :
                             h2->Fill(CFDt);;
                             break;
                         case 3 :
                             h3->Fill(CFDt);;
                             break;
                         case 4 :
                             h4->Fill(CFDt);;
                             break;
                         case 5 :
                             h5->Fill(CFDt);;
                             break;
                         case 6 :
                             h6->Fill(CFDt);;
                             break;
                         case 7 :
                             h7->Fill(CFDt);;
                             break;
                         case 8 :
                             h8->Fill(CFDt);;
                             break;
                         case 9 :
                             h9->Fill(CFDt);;
                             break;
                         default:
                             cout << " - - CFD Error!"<<endl;
                     }

                  }
             //cout<<"--Waveform NUM: "<<j<<"	peak value:"<<pv-v_baseline<<"mV	%5-value:"<<v1-v_baseline<<"mV	%95-value:"<<v2-v_baseline<<"mV	v-cfd:"<<v_cfd-v_baseline<<"mV	baseline:"<<v_baseline<<"mV	point num:"<<point_num[0]<<"	tt:"<<CFDt<<"ns"<<endl;
             //cout<<"--Waveform NUM: "<<j<<"	peak value:"<<peak1[0]-v_baseline1<<"mV	CFDt:"<<CFDt<<"ns"<<endl;
         }

    }
        TCanvas *c2 = new TCanvas("C3*3","fit",200,10,2000,500);
        c2->Divide (5,2);
        //c2->SetFillColor(42);
        c2->SetGrid();

        TF1 *g1 = new TF1("g1","gaus",0,5);
        double par[27];
        c2->cd(1);
        h1->Draw();
        h1->Fit(g1); g1->GetParameters(&par[0]);
        //cout<<"----------------"<<par[2]<<endl;
        //h1->SetOptFit("1");
        h1->SetLineWidth(2);
        h1->SetLineColor(1);
        h1->GetXaxis()->SetTitle("CFDt [ns]");
        h1->GetYaxis()->SetTitle("Count");
        h1->GetYaxis()->SetTitleOffset(1.3);
        gStyle->SetOptFit(1011);
        c2->cd(2); h2->Draw();h2->Fit(g1); g1->GetParameters(&par[3]);
        c2->cd(3); h3->Draw();h3->Fit(g1); g1->GetParameters(&par[6]);
        c2->cd(4); h4->Draw();h4->Fit(g1); g1->GetParameters(&par[9]);
        c2->cd(5); h5->Draw();h5->Fit(g1); g1->GetParameters(&par[12]);
        c2->cd(6); h6->Draw();h6->Fit(g1); g1->GetParameters(&par[15]);
        c2->cd(7); h7->Draw();h7->Fit(g1); g1->GetParameters(&par[18]);
        c2->cd(8); h8->Draw();h8->Fit(g1); g1->GetParameters(&par[21]);
        c2->cd(9); h9->Draw();h9->Fit(g1); g1->GetParameters(&par[24]);

        double CFDFraction[9]={0.1,.2,.3,.4,.5,.6,.7,.8,.9},\
               parF[9]={1000*par[2],1000*par[5],1000*par[8],1000*par[11],1000*par[14],1000*par[17],1000*par[20],1000*par[23],1000*par[26]};
        TGraph *gF=new TGraph(9,CFDFraction,parF);
        c2->cd(10);
        gF->Draw("AC*");
        gF->SetTitle("CFD Fraction ");
        gF->GetXaxis()->SetTitle("CFD Fraction");
        gF->GetYaxis()->SetTitle("Resolution [ps]"); //

}


void T_CFD2(){

        TString file1,file2;
        double CFD = 0.1;                    //CFD Fraction [0-1]
        int FillCounter =0;
        int j;
        double *peak1,*peak2,PeakV1,PeakV2,PeakP1,PeakP2;
        double CFDt1,CFDt2,CFDt,peak_fit1[2],peak_fit2[2],v_baseline1,v_baseline2,q;
        static double waveform_point1[STEP],waveform_point2[STEP];
        TH1F *h1=new TH1F("Cap10_0.1","CFD-0.1",500,-1,1);//0-1ns
        TH1F *h2=new TH1F("Cap10_0.2","CFD-0.2",500,-1,1);//0-1ns
        TH1F *h3=new TH1F("Cap10_0.3","CFD-0.3",500,-1,1);//0-1ns
        TH1F *h4=new TH1F("Cap10_0.4","CFD-0.4",500,-1,1);//0-1ns
        TH1F *h5=new TH1F("Cap10_0.5","CFD-0.5",500,-1,1);//0-1ns
        TH1F *h6=new TH1F("Cap10_0.6","CFD-0.6",500,-1,1);//0-1ns
        TH1F *h7=new TH1F("Cap10_0.7","CFD-0.7",500,-1,1);//0-1ns
        TH1F *h8=new TH1F("Cap10_0.8","CFD-0.8",500,-1,1);//0-1ns
        TH1F *h9=new TH1F("Cap10_0.9","CFD-0.9",500,-1,1);//0-1ns
        for(j=0;j<N;j++){
             char Jstr[5]={0};
             sprintf(Jstr,"%05d",j);
             file1= title+Jstr;//to_string(Jstr);    //  BV170_130V_19999Wfm_Ch1
             file2= title2+Jstr;//to_string(Jstr);
             cout<<" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  B e g i n  "<<j+1<<endl;
             get_waveform_data(file1,waveform_point1,STEP);
             get_waveform_data(file2,waveform_point2,STEP);
             //cout<<" - - - - - - - - - - - - - - - - - - - - - - - - - check-1"<<endl;
             peak1 = get_peak_parameters(waveform_point1,peak_fit1,THRESHOLD,STEP);
             PeakV1 = peak1[0];  PeakP1 = peak1[1];
             peak2 = get_peak_parameters(waveform_point2,peak_fit2,THRESHOLD,STEP);
             PeakV2 = peak2[0];  PeakP2 = peak2[1];
             cout<<"--PeakV1 "<<PeakV1<<"mV	 PeakV2 "<<PeakV2<<"mV"<<endl;
             //cout<<"--PeakFit  :"<<peak_fit1[0]<<" mV	PeakV2:"<<peak_fit2[0]<<" mV"<<endl;
             v_baseline1 = get_baseline(waveform_point1,BASELINE_START,BASELINE_STOP);
             v_baseline2 = get_baseline(waveform_point2,BASELINE_START,BASELINE_STOP);

             //if(PeakV1>THRESHOLD && PeakV1<700 && PeakP1*DELTA_TIME>6 && PeakP1*DELTA_TIME<8 && PeakV2>THRESHOLD && PeakV2<700){
             if(PeakV1>THRESHOLD && PeakV1<700 && PeakV2>THRESHOLD && PeakV2<700){
             //if(PeakV1>THRESHOLD &&  PeakV2>THRESHOLD ){
                   FillCounter++;
                   cout<<"---------------------------------------------------------------------- Fill counter "<<FillCounter<<endl;
                  for(CFD = 0.1;CFD<0.91;CFD+=0.1){
                     //CFD = 0.8;
                     CFDt1 = get_CFDt2(waveform_point1,v_baseline1,CFD,PeakV1,PeakP1);
                     CFDt2 = get_CFDt2(waveform_point2,v_baseline2,CFD,PeakV2,PeakP2);
                     CFDt = CFDt1-CFDt2;

                     int intCFD = CFD*10+0.1;
                     cout<<" - - - intCFD "<<intCFD<<"  CFD "<<CFD<<endl;

                     switch (intCFD)
                     {
                         case 1 :
                             h1->Fill(CFDt);;
                             break;
                         case 2 :
                             h2->Fill(CFDt);;
                             break;
                         case 3 :
                             h3->Fill(CFDt);;
                             break;
                         case 4 :
                             h4->Fill(CFDt);;
                             break;
                         case 5 :
                             h5->Fill(CFDt);;
                             break;
                         case 6 :
                             h6->Fill(CFDt);;
                             break;
                         case 7 :
                             h7->Fill(CFDt);;
                             break;
                         case 8 :
                             h8->Fill(CFDt);;
                             break;
                         case 9 :
                             h9->Fill(CFDt);;
                             break;
                         default:
                             cout << " - - CFD Error!"<<endl;
                     }

                  }
             //cout<<"--Waveform NUM: "<<j<<"	peak value:"<<pv-v_baseline<<"mV	%5-value:"<<v1-v_baseline<<"mV	%95-value:"<<v2-v_baseline<<"mV	v-cfd:"<<v_cfd-v_baseline<<"mV	baseline:"<<v_baseline<<"mV	point num:"<<point_num[0]<<"	tt:"<<CFDt<<"ns"<<endl;
             //cout<<"--Waveform NUM: "<<j<<"	peak value:"<<peak1[0]-v_baseline1<<"mV	CFDt:"<<CFDt<<"ns"<<endl;
         }

    }
        TCanvas *c2 = new TCanvas("C3*3","fit",200,10,2000,500);
        c2->Divide (5,2);
        //c2->SetFillColor(42);
        c2->SetGrid();

        TF1 *g1 = new TF1("g1","gaus",0,5);
        double par[27];
        c2->cd(1);
        h1->Draw();
        h1->Fit(g1); g1->GetParameters(&par[0]);
        //cout<<"----------------"<<par[2]<<endl;
        //h1->SetOptFit("1");
        h1->SetLineWidth(2);
        h1->SetLineColor(1);
        h1->GetXaxis()->SetTitle("CFDt [ns]");
        h1->GetYaxis()->SetTitle("Count");
        h1->GetYaxis()->SetTitleOffset(1.3);
        c2->cd(2); h2->Draw();h2->Fit(g1,"Q"); g1->GetParameters(&par[3]);
        c2->cd(3); h3->Draw();h3->Fit(g1,"Q"); g1->GetParameters(&par[6]);
        c2->cd(4); h4->Draw();h4->Fit(g1,"Q"); g1->GetParameters(&par[9]);
        c2->cd(5); h5->Draw();h5->Fit(g1,"Q"); g1->GetParameters(&par[12]);
        c2->cd(6); h6->Draw();h6->Fit(g1,"Q"); g1->GetParameters(&par[15]);
        c2->cd(7); h7->Draw();h7->Fit(g1,"Q"); g1->GetParameters(&par[18]);
        c2->cd(8); h8->Draw();h8->Fit(g1,"Q"); g1->GetParameters(&par[21]);
        c2->cd(9); h9->Draw();h9->Fit(g1,"Q"); g1->GetParameters(&par[24]);gStyle->SetOptFit(1011);

        double CFDFraction[9]={0.1,.2,.3,.4,.5,.6,.7,.8,.9},\
               parF[9]={1000*par[2],1000*par[5],1000*par[8],1000*par[11],1000*par[14],1000*par[17],1000*par[20],1000*par[23],1000*par[26]};
        TGraph *gF=new TGraph(9,CFDFraction,parF);
        c2->cd(10);
        gF->Draw("AC*");
        gF->SetTitle("CFD Fraction ");
        gF->GetXaxis()->SetTitle("CFD Fraction");
        gF->GetYaxis()->SetTitle("Resolution [ps]"); //

        string CName = title;
        cout<<"qian: "<<title<<endl;
        replace(CName.begin(),CName.end(),'/','-');
        cout<<"hou: "<<CName<<endl;
        c2->SaveAs(Form("CFD_%s.root",CName.data()));




}


double get_CFDt1(double *waveform,double v_baseline,double CFD,double PeakV,double PeakP){    //CFD, not fit waveform. For WF2. Add Noise For No Noise Signal.

        double CFDt;
        double v_cfd = CFD*(PeakV-v_baseline)+v_baseline;
        int i_cfd;
        //cout<<" V_CFD "<<v_cfd<<endl;
        double Noise = 1.25;  //noise [mV]
        double waveformNoise = gRandom->Gaus()*Noise;
        cout<<"waveformNoise    "<<waveformNoise<<endl;
        double slope;

        for(int i = PeakP;i>0;i--){
            if (waveform[i]+waveformNoise < v_cfd){
                i_cfd = i;
                CFDt = i*DELTA_TIME;
                slope = abs(waveform[i-1] - waveform[i+1])/(2*DELTA_TIME) ;
                break;
            }
        }

        double NoiseJitter = 0;
        float Jit =  0;
        if (slope !=0 )
          {
            NoiseJitter = Noise/slope; // in ns
            Jit  = gRandom->Gaus(0, NoiseJitter);
          }
        CFDt = CFDt+Jit;
        cout<<"Noise--Jit    "<<Jit<<endl;

        //cout<<"-----------------CFDt "<<CFDt<<endl;

    return(CFDt);
}


double get_CFDt2(double *waveform,double v_baseline,double CFD,double PeakV,double PeakP){    //CFD,fit waveform

        double CFDt;
        double v_cfd = CFD*(PeakV-v_baseline)+v_baseline;
        int i_cfd, fitPointNum = 6;   //CFD fit range (i_cfd-fitPointNum/2+1   i_cfd+fitPointNum/2)
        double fitT1,fitT2,fitStatus=0;
        //cout<<" V_CFD "<<v_cfd<<endl;

        vector <double> value;              //
        vector <double> t;

        for(int i = PeakP;i>0;i--){
            if (waveform[i] < v_cfd){ i_cfd = i; break; }
        }

        for(int j = i_cfd-fitPointNum/2;j<i_cfd+fitPointNum/2+1;j++){
            value.push_back(waveform[j]);
            t.push_back(j*DELTA_TIME);
        }

        fitT1 = (i_cfd-fitPointNum/2    )*DELTA_TIME;
        fitT2 = (i_cfd+fitPointNum/2+1  )*DELTA_TIME;
        cout<<"--Fit range [ns]  "<<fitT1<<" - "<<fitT2<<endl;

        //TCanvas *cCFDt = new TCanvas("CCFDt","fitCFDt",200,10,2000,500);
        TGraph *g=new TGraph(fitPointNum,&t[0],&value[0]);
        TF1 *f=new TF1("f","pol3",fitT1,fitT2);
        //TF1 *f = new TF1("f", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", fitT1, fitT2);
        //TF1 *f = new TF1("f", "[0]+[1]*x+[2]*x*x", fitT1, fitT2);
        //TF1 *f = new TF1("f", "[0]+[1]*x", fitT1, fitT2);
        g->Fit(f,"QR+");
        if(f->GetChisquare()>1e-4) fitStatus=1;
        CFDt = f->GetX(v_cfd, fitT1, fitT2, 1e-14);
        //cout<<"-----------------CFDt "<<CFDt<<endl;
        //g->Draw("AC*");
    delete g;
    delete f;
    return(CFDt);
}


/*
void Baseline_distribution(){
        TString file0,num;
        int i;
        double pp,peak_par[2],v_b;
        TH1F *h=new TH1F("h","Baseline Distribution",400,-2,2);
        for(i=0;i<N;i++){
                num=itoa(i,5);
                file0=title+"_"+num;
                get_peak_parameters(file0,peak_par,THRESHOLD,STEP);
                pp=peak_par[1];
                v_b=get_baseline(file0,BASELINE_START,BASELINE_STOP);
                if(pp>0){
                        h->Fill(-1*v_b);
                        cout<<"NUM: "<<i<<" Baseline:"<<-1*v_b<<"mV"<<endl;
                }
        }
        TCanvas *c5 = new TCanvas("c5","fit",200,10,1000,500);
        //c5->SetFillColor(42);
        c5->SetGrid();
        gStyle->SetOptFit(1100);
        //gPad->SetLogy(1);
        h->Draw();
        h->SetLineWidth(2);
        h->SetLineColor(1);
        h->GetXaxis()->SetTitle("U/mV");
        h->GetYaxis()->SetTitle("count");
        h->GetYaxis()->SetTitleOffset(1.3);
}

*/

void Charge_distribution(double gate_width){                 //charge spectrum

        TH1F *h=new TH1F("h","charge spectrum",100,0,20000);
        string file1;
        int j;
        double *peak1,PeakV1,PeakP1,peak_fit1[2],Q,v_baseline1;
        static double waveform_point1[STEP];

        for(j=0;j<N;j++){
            cout<<" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  B e g i n  "<<j+1<<endl;
            char Jstr[5]={0};
            sprintf(Jstr,"%05d",j);
            int Vt = THRESHOLD;                            //mV
            file1= title2+Jstr;//to_string(Jstr);    //  BV170_130V_19999Wfm_Ch1
            //file2= title2+Jstr;//to_string(Jstr);
            //file1=title+to_string(j+1)+"Wfm_Ch2";    //

            get_waveform_data(file1,waveform_point1,STEP);

            peak1 = get_peak_parameters(waveform_point1,peak_fit1,THRESHOLD,STEP);
            PeakV1 = peak1[0];  PeakP1 = peak1[1];
            v_baseline1 = get_baseline(waveform_point1,BASELINE_START,BASELINE_STOP);

            int gateL = PeakP1-gate_width/DELTA_TIME/2;
            int gateR = PeakP1+gate_width/DELTA_TIME/2;
            //cout<<"gateL "<<gateL<<" gateR "<<gateR<<" i"<<endl;

            int count_point = 0;
            for(int i=PeakP1;i<gateR+1 ;i++){
                Q+=waveform_point1[i];
                count_point++;
                if(waveform_point1[i]<v_baseline1) {cout<<"gateR:  "<<i*DELTA_TIME<<" ns"<<endl; break;}
            }

            for(int i=PeakP1-1;i>gateL ;i--){
                Q+=waveform_point1[i];
                count_point++;
                if(waveform_point1[i]<v_baseline1) {cout<<"gateL:  "<<i*DELTA_TIME<<" ns"<<endl; break;}
            }

            Q-=count_point*v_baseline1;
            Q*=20*DELTA_TIME;   // 1000/50=20  mA->uA     Q[fC]

            if(Q>Q_THRESHOLD && PeakV1>Vt){// && PeakV1< 355){
                 cout<<"N: "<<j<<" -- Charge:  "<<Q<<" fC"<<endl;
                 h->Fill(Q);
             }
        }
        TCanvas *cC = new TCanvas("cC","Charge_distribution",200,10,1000,500);
        //c6->SetFillColor(42);
        cC->SetGrid();
        gStyle->SetOptFit(1100);
        //gPad->SetLogy(1);
        h->Draw();
        h->SetLineWidth(2);
        h->SetLineColor(1);
        h->GetXaxis()->SetTitle("Charge [fC]");
        h->GetYaxis()->SetTitle("count");
        h->GetYaxis()->SetTitleOffset(1.3);

        string CName = file1;
        cout<<"qian: "<<title<<endl;
        replace(CName.begin(),CName.end(),'/','-');
        cout<<"hou: "<<CName<<endl;
        TFile *fout = new TFile(Form("Charge_%s.root",CName.data()),"recreate" );
        //TFile *fout = new TFile(Form("Charge.root",CName.data()) );
        h->Write("h");
}

void Peak_distribution(){                 //peak spectrum


        TH1F *h=new TH1F("h","peak spectrum",100,0,500);
        string file1;
        int j;
        double *peak1, peak_fit1[2];
        static double waveform_point1[STEP];

        for(j=0;j<N;j++){
            cout<<" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  B e g i n  "<<j+1<<endl;
            char Jstr[5]={0};
            sprintf(Jstr,"%05d",j);
            file1= title+Jstr;//to_string(Jstr);    //  BV170_130V_19999Wfm_Ch1
            //file2= title2+Jstr;//to_string(Jstr);
            //file1=title+to_string(j+1)+"Wfm_Ch2";    //

            get_waveform_data(file1,waveform_point1,STEP);

            peak1 = get_peak_parameters(waveform_point1,peak_fit1,THRESHOLD,STEP);

            if(peak1[0]>THRESHOLD){
            //if(peak1[0]>THRESHOLD && peak1[1]*DELTA_TIME>6 && peak1[1]*DELTA_TIME<9){
            //if(peak1[0]>THRESHOLD && peak1[1]*DELTA_TIME>9 && peak1[1]*DELTA_TIME<13){
                 cout<<"N: "<<j<<" -- Peak:  "<<peak1[0]<<" mV"<<endl;
                 h->Fill(peak1[0]);
             }
        }
        TCanvas *c6 = new TCanvas("c6","Peak_distribution",200,10,1000,500);
        //c6->SetFillColor(42);
        c6->SetGrid();
        gStyle->SetOptFit(1100);
        //gPad->SetLogy(1);
        h->Draw();
        h->SetLineWidth(2);
        h->SetLineColor(1);
        h->GetXaxis()->SetTitle("U [mV]");
        h->GetYaxis()->SetTitle("count");
        h->GetYaxis()->SetTitleOffset(1.3);

        string CName = file1;
        cout<<"qian: "<<title<<endl;
        replace(CName.begin(),CName.end(),'/','-');
        cout<<"hou: "<<CName<<endl;
        TFile *fout = new TFile(Form("Peak_%s.root",CName.data()),"recreate" );
        //TFile *fout = new TFile(Form("Peak.root",CName.data()) );
        h->Write("h");
}


void Noise_distribution(double win_start,double win_stop){

        TH1F *hN=new TH1F("Noise","Noise",50,-10,10);//-10-10mV

        TString file1;
        static double waveform_point1[STEP];
        int j,b,s;

        for(j=0;j<N;j++){
            cout<<" * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  B e g i n  "<<j+1<<endl;
            file1=title+to_string(j+10)+"Wfm_Ch2";    //
            get_waveform_data(file1,waveform_point1,STEP);

            b=(int)(win_start/DELTA_TIME);
            s=(int)(win_stop/DELTA_TIME);

            for (int i=b;i<s;i++){
                hN->Fill(waveform_point1[i]);
            }
        }

        TCanvas *cN = new TCanvas("CN","Noise",200,10,600,400);
        TF1 *g1 = new TF1("g1","gaus",-10,10);
        double par[3];
        hN->Draw();
        hN->Fit(g1); g1->GetParameters(&par[0]); gStyle->SetOptFit(1011);
        hN->SetLineWidth(2);
        hN->SetLineColor(1);
        hN->GetXaxis()->SetTitle("Noise [mV]");
        hN->GetYaxis()->SetTitle("Count");
        hN->GetYaxis()->SetTitleOffset(1.3);

}
//*/

double *get_peak_parameters(double *waveform_point,double *peak_fit,int threshold,int step){

    int i,j=0;
    const int n = 20;
    double temp2,x[n]={0},y[n]={0};
    static double peak[2];
    peak[0] = 0;
    peak[1] = 0;

    for(peak[0] = 0,i=0;i<step;i++){
        temp2 = waveform_point[i];
        //cout<<min<<endl;
        if(temp2 > peak[0]){
           peak[0] = temp2;
           peak[1]=i;
        }
    }
    //cout<<"--Peak position-max:"<<peak[1]*DELTA_TIME<<"ns	peak-max:"<<peak[0]<<"mV"<<endl;
    /*
    if(peak[0]>threshold){
        //cout<<"peak_position = "<<peak[1]<<" ns"<<endl;
        for(i=peak[1]-n/2;i<peak[1]+n/2;i++){
            if(i>=step){
                    break;
            }
            if(i<0){
                    continue;
            }
            x[j]=i;
            y[j]=waveform_point[i];
            //cout<<x[j]<<"   "<<y[j]<<endl;
            j++;
        }
        //TCanvas *c2 = new TCanvas("Peakfit","Peakfit",200,10,1000,500);
        //c2->SetGrid();
        TGraph *g=new TGraph(n,x,y);
        TF1 *f=new TF1("f","gaus",peak[1]-10,peak[1]+11);
        //g->Draw("AC*");
        g->Fit(f,"R+");
        f->GetParameters(&peak_fit[0]);
        peak_fit[1]*=1.0*DELTA_TIME;
        peak_fit[0]*=1;
        cout<<"--Peak position-fit:"<<peak_fit[1]<<"ns  peak-fit:"<<peak_fit[0]<<"mV"<<endl;
            delete g;
            delete f;
    }
    else{
        peak_fit[0]=peak_fit[1]=0;
    }  //  */

    return(peak);
}


double get_baseline(double *waveform_point,double win_start,double win_stop){

        int b,s,width=0;
	b=(int)(win_start/DELTA_TIME);
	s=(int)(win_stop/DELTA_TIME);

	double sum=0.0,baseline;

        for (int i=b;i<s;i++){
            sum+=waveform_point[i];
            width++;
	}

        baseline=sum/width;
	return(baseline);
}


void get_waveform_data(TString file,double *waveform,int step){    //for .csv file
        ifstream in;
        cout<<file<<endl;
        in.open(file+".csv");
        int i;
        double t[STEP];

        //line 1 to line 5 are comments
        char tmp[255];
        for( int lineNo=0; lineNo<5; lineNo++)
        {
          in.getline(tmp, 255);
          //cout<<"tmp  -----  "<<tmp<<endl;
        }// */

        string line;
        for(i=0;i<1000;i++){
            getline(in,line);//读取每行数据
            string number;
            istringstream readstr(line); //string数据流化
            //将一行数据按'，'分割
            for(int j = 0;j < 2;j++){ //可根据数据的实际情况取循环获取
                getline(readstr,number,','); //循环读取数据
                if (j == 1) waveform[i] = (atof(number.c_str())); // str to float
            }
            waveform[i] *= -1000;//
            //cout<<"-------------waveform  "<<i<<"   "<<waveform[i]<<endl;
            t[i]=i*DELTA_TIME;
        }

        //if(file == "BV170_130V_1000Wfm_Ch1"){
            TGraph *g=new TGraph(step,t,waveform);
            TCanvas *c = new TCanvas("Waveform","Waveform",200,10,1000,500);
            c->SetGrid();
            g->Draw("AC*");
            g->SetTitle("Waveform of "+file);
            g->GetXaxis()->SetTitle("t/ns");
            g->GetYaxis()->SetTitle("V/mV");
            //g->GetYaxis()->SetTitle("I/uA");
        //}

        in.close();
}
//*/

/*
void get_waveform_data(TString file,double *waveform,int step){  // for .txt file
	ifstream in;
        cout<<file<<endl;
        in.open(file+".txt");
	int i;
        double t[STEP];
        double temp1;

        // line 1 to line 4 are comments
        char tmp[255];
        for( int lineNo=0; lineNo<4; lineNo++)
        {
          in.getline(tmp, 255);
        }

        for(i=0;i<step;i++){
                in>>temp1>>temp1>>temp1>>temp1>>temp1>>temp1>>waveform[i]>>temp1>>temp1; //Signal after BB amplifier
                //in>>temp1>>waveform[i]>>temp1>>temp1>>temp1>>temp1>>temp1>>temp1>>temp1;  //Sensor current
                waveform[i] *= -1;//
                //cout<<"-------------waveform  "<<i<<"   "<<waveform[i]<<endl;
                t[i]=i*DELTA_TIME;
	}

        //if(file == "BV170_130V_2Wfm_Ch1"){
            TGraph *g=new TGraph(step,t,waveform);
            TCanvas *c = new TCanvas("Waveform","Waveform",200,10,1000,500);
            c->SetGrid();
            g->Draw("AC*");
            g->SetTitle("Waveform of "+file);
            g->GetXaxis()->SetTitle("t/ns");
            g->GetYaxis()->SetTitle("V/mV");
            //g->GetYaxis()->SetTitle("I/uA");
        //}

        in.close();
}
//*/


//





