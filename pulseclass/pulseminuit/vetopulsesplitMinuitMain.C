#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TRint.h"
#include "TProof.h"
#include "TCut.h"
#include "TLegend.h"
#include "TMinuit.h"

using namespace std;
#include "vetopulsesplit.hh"
#include "./odselector/odselector.h"
#define StartFit

TRint* theApp;
string Time = "Mar22PM";
string dirname="/darkside/users/hqian/pulse_splitter/";

void Process(TChain* chain, TString label){
  TString option = label;
  chain->SetProof();
  TProof* pr = TProof::Open("workers=2");
  pr->SetParameter("PROOF_Packetizer","TPacketizer");
  pr->SetParameter("PROOF_MaxSlavesPerNode",8);
  chain->Process("./odselector/odselector.C+",option.Data(),-1,0);
}

#ifndef __CINT__
int main(int argc, char** argv){
  theApp = new TRint("App",&argc,argv,NULL,0);
  int start, end,tag;
  if(theApp->Argc() == 1)
    {
      tag = 1;
    }
  else if(theApp->Argc() == 3)
    {
      tag=3;
      start = atoi(theApp->Argv(1));
      end = atoi(theApp->Argv(2));
    }
  else{
    cout<<"Usage: ./reconMain startfile endfile "<<endl;
    cout<<"Usage: ./reconMain "<<endl;
    return 0;
  }

  vetopulsesplit *split = new vetopulsesplit;
  split->SetMCinputdir(dirname+"MCData/");
  split->SetRealinputdir(dirname+"VetoData/");
  split->Setoutputdir(dirname);  
  if(tag==3)
    {
      //      split->ODTree_DataAnalysis(start,end);      
      TChain *odtree =new TChain("odtree");
      split->Readdatafile(odtree,start,end);
      Int_t nEntries = odtree->GetEntries();
      cout<<"nEntries= "<<nEntries<<endl;
      string label = dirname+"VetoData/";
      Process(odtree,label);
    }
  else if(tag==1){
    split->Recon_DataAnalysis();
    //    string RealDataFile = split->GetRealinputdir() +"PulseSplitRealEnergy.root";
    string RealDataFile = split->GetRealinputdir() +"PulseSplitRealEnergy_MoreBins.root";
    cout<<"Fitting "<<RealDataFile<<endl;
    TFile *RealFile = new TFile(RealDataFile.c_str());    
    TH1F *FullSpectrum = (TH1F*) RealFile->Get("FullSpectrum");
    TLegend *leg = new TLegend(0.3,0.55,0.7,0.9);
    TCanvas *c1 = new TCanvas("c1", "Full Spectrum with Fit", 1200, 600);
    c1->SetLogy();
    c1->cd();
    FullSpectrum->Draw();
  
    double startrange=25;
    double endrange = 1600;
    int startbin = FullSpectrum->FindBin(startrange);
    int endbin = FullSpectrum->FindBin(endrange);
    Double_t integralsum = FullSpectrum->Integral(startbin,endbin);
    const Int_t parnums =15;
    const Int_t isotopes = 6;
    string Source[] = {"C14","C14_1st","C14_2nd","Co60","Co57","K40","Th232","U238"};
    Double_t isotope_activity[isotopes] = {223268845,223268845,223268845,99473576,24773624,47584434};
    Double_t isotope_fraction[isotopes];
    split->Isotope_Fraction(isotopes,isotope_activity,isotope_fraction);
    Double_t integral_sum[isotopes];
    for(int i=0; i<isotopes; i++)
      integral_sum[i] = isotope_fraction[i]*integralsum;
    
    Double_t new_par[parnums];
    
    Double_t binw             = FullSpectrum->GetBinWidth(1);
    Double_t ly_mean          = 0.55;
    Double_t spe_var          = 0.14;
    Double_t ly_var           = 0.00642;
    Double_t constant         = 0.;
    Double_t basel_mean       = 0.386;
    Double_t basel_var        = 1.;
    Double_t threshold        = 1;
    Double_t window           = 2e-7;

    Double_t parvaluesamples[]={binw,ly_mean,spe_var,ly_var,constant,basel_mean,basel_var,threshold,window,
				//integralsum,integralsum,integralsum,integralsum,integralsum,integralsum,integralsum};
				integral_sum[0],integral_sum[1],integral_sum[2],integral_sum[3],integral_sum[4],integral_sum[5]};
    string parnamesamples[]={"Bin Width","LY Mean [PE/keV]","Rel SPE Var","Rel LY Var",
			     "Constant","Baseline Mean [p.e]","Baseline Var","Threshold","Window Width [s]",
			     "C14 Amplitude","Co60 Amplitude","Co57 Amplitude","K40 Amplitude","Th232 Amplitude","U238 Amplitude"};
    Double_t parlowlimits[parnums] = {0,0,0,-100,-100,-100,-100,-100,0,0,0,0,0,0,0};
    Double_t paruplimits[]  = {10,1.,0.14,100,100,100,100,100,1e-6,integralsum,integralsum,integralsum,integralsum,integralsum,integralsum};     //1,1,1,1,1,1};
    Double_t step[parnums] = {1,0.01,0.01,0.0001,0.0001,0.001,0.1,0.4,2.e-9,1,1,1,1,1,1};

    TMinuit *gMinuit = new TMinuit(parnums);
    gMinuit->SetFCN(fcn);
    Double_t arglist[];
    Int_t ierflag = 0;
    arglist[0]=1;
    gMinuit->mnexcm("SET ERR",arglist,1,ierflag);
    for(int i=0, i<parnums; i++)
      gMinuit->mnparm(i,parnamesamples[i],parvaluesamples[i],step[i],parlowlimits[i],paruplimits[i],ierflag);	     
    arglist[0]=500;
    arglist[1]=0;
    gMinuit->mnexcm("MIGRAD",arglist,2,ierflag);
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

    TF1* Fit_C14 = new TF1("Fit_C14",split,&vetopulsesplit::C14Fit, startrange, endrange, parnums,"vetopulsesplit","C14Fit");        
    TF1* Fit_C14_1st = new TF1("Fit_C14_1st",split,&vetopulsesplit::C14Fit_1st , startrange, endrange, parnums,"vetopulsesplit","C14Fit_1st"); 
    TF1* Fit_C14_2nd = new TF1("Fit_C14_2nd", split,&vetopulsesplit::C14Fit_2nd, startrange, endrange, parnums,"vetopulsesplit","C14Fit_2nd"); 
    //   FitC14 = new TF1("FitC14",C14Fit,start[0],end[0],parnums); //total C14 fit function
    TF1* Fit_Co60  = new TF1("Fit_Co60",split,&vetopulsesplit::Co60Fit,startrange,endrange,parnums,"vetopulsesplit","Co60Fit");
    TF1* Fit_Co57  = new TF1("Fit_Co57",split,&vetopulsesplit::Co57Fit,startrange,endrange,parnums,"vetopulsesplit","Co57Fit");
    TF1* Fit_K40   = new TF1("Fit_K40",split,&vetopulsesplit::K40Fit,startrange,endrange,parnums,"vetopulsesplit","K40Fit");
    //  TF1* Fit_Tl208 = new TF1("Fit_Tl208",split,&vetopulsesplit::Tl208Fit,startrange,endrange,parnums,"vetopulsesplit","Tl208Fit");
    TF1* Fit_Th232 = new TF1("Fit_Th232",split,&vetopulsesplit::Th232Fit,startrange,endrange,parnums,"vetopulsesplit","Th232Fit");
    TF1* Fit_U238  = new TF1("Fit_U238",split,&vetopulsesplit::U238Fit,startrange,endrange,parnums,"vetopulsesplit","U238Fit");
    TF1* Fit_Total  = new TF1("Fit_Total",split,&vetopulsesplit::TotalFit,startrange,endrange,parnums,"vetopulsesplit","TotalFit");
    TObjArray Fitlist(0);
    Fitlist.Add(Fit_C14);
    Fitlist.Add(Fit_C14_1st);
    Fitlist.Add(Fit_C14_2nd);
    Fitlist.Add(Fit_Co60);
    Fitlist.Add(Fit_Co57);
    Fitlist.Add(Fit_K40);
    // Fitlist.Add(Fit_Tl208);
    Fitlist.Add(Fit_Th232);
    Fitlist.Add(Fit_U238);
    //  Fitlist.Add(Fit_Total);
    vector<int> linecolor = split->Colors();
    for(int i=0; i<parnums; i++)
      {
	Fit_Total->SetParName(i,parnamesamples[i].c_str());
	if(i==4 || i == 0 || i==2) Fit_Total->FixParameter(i,parvaluesamples[i]);
	else{
	  Fit_Total->SetParameter(i, parvaluesamples[i]);
	  Fit_Total->SetParLimits(i,parlowlimits[i],paruplimits[i]);
	}
      }
    Fit_Total->SetNpx(2000);
    //  FitTotal->SetLineColor(kRed);
    Fit_Total->SetLineColor(linecolor.at(Fitlist.GetEntries()+1));  
#ifdef StartFit
    FullSpectrum->Fit(Fit_Total,"RVEM");
#endif
    Fit_Total->GetParameters(new_par);

    for(Int_t i=0; i<Fitlist.GetEntries(); i++)
      {
	TF1* tempfit = dynamic_cast<TF1*>(Fitlist.At(i));
	cout<<tempfit->GetName()<<endl;
	tempfit->SetLineColor(linecolor.at(i));
	tempfit->SetParameters(new_par);
	tempfit->Draw("SAME");
	//	leg->AddEntry(tempfit,Source[i].c_str(),"l");
	leg->AddEntry(tempfit,tempfit->GetName(),"l");
	cout<<"4 "<<tempfit->GetName()<<"  "<<Fitlist.GetSize()<<" "<<Fitlist.GetEntries()<<endl;
      }
#ifndef StartFit
    Fit_Total->Draw("same");
#endif
    leg->AddEntry(Fit_Total,"Total Fit","l");
    leg->Draw();

    Double_t chi2 = Fit_Total->GetChisquare();
    Int_t ndf = Fit_Total->GetNDF();
    cout<<"chi2/ndf= "<<chi2/(1.0*ndf)<<"\t Probability="<<Fit_Total->GetProb()<<endl;
    string outdir=split->Getoutputdir();
    string output = outdir + "Pulse_Splitter_"+Time+".root";
    TFile outfile(output.c_str(), "RECREATE");
    split->GetHistArray().Write();
    Fitlist.Write();
    Fit_Total->Write();
    FullSpectrum->Write();
    c1->Write();
    outfile.Write();
    outfile.Close();
  
    delete split;
    cout<<"Successfully finished the app."<<endl;
  }
  return 1;  
}
#endif /*__CINT__*/




