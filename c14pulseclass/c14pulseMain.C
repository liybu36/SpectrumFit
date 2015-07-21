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
#include "TFitResultPtr.h"

using namespace std;
#include "c14pulse.hh"

TRint* theApp;
string Time = "Jun24";
string dirname="/darkside/users/hqian/pulse_splitter/";

#ifndef __CINT__
int main(int argc, char** argv){
  theApp = new TRint("App",&argc,argv,NULL,0);
  int start, end,tag;
  if(theApp->Argc() == 1)
    {
      tag = 1;
      //  start = atoi(theApp->Argv(1));
      //  end = start;
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

  c14pulse *split = new c14pulse;
  split->SetMCinputdir(dirname+"SDMCData/");
  split->SetRealinputdir(dirname+"VetoData/");
  split->Setoutputdir(dirname);  
  /*  if(tag==3)
    {
#include "./odselector/odselector.h"
      void Process(TChain* chain, TString label){
	TString option = label;
	chain->SetProof();
	TProof* pr = TProof::Open("workers=2");
	pr->SetParameter("PROOF_Packetizer","TPacketizer");
	pr->SetParameter("PROOF_MaxSlavesPerNode",8);
	chain->Process("./odselector/odselector.C+",option.Data(),-1,0);
      }
      //      split->ODTree_DataAnalysis(start,end);      
      TChain *odtree =new TChain("odtree");
      split->Readdatafile(odtree,start,end);
      Int_t nEntries = odtree->GetEntries();
      cout<<"nEntries= "<<nEntries<<endl;
      string label = dirname+"VetoData/";
      Process(odtree,label);
*/      /*
      TChain *odtree =new TChain("odtree");
      split->Readdatafile(odtree,start,end);
      Int_t nEntries = odtree->GetEntries();
      cout<<"nEntries= "<<nEntries<<endl;      
      double startrange=0;
      double endrange = 2100;
      TH1F* FullSpectrum = new TH1F("FullSpectrum","lsv cluster fixed width charge full spectrum;charge [PE];Counts",1000,startrange,endrange);
      TCut TimingCut = "lsv_cluster_start_ns.fArray>3770 && lsv_cluster_start_ns.fArray<3786";
      TCut MultiCut = "lsv_cluster_height.fArray/lsv_cluster_max_multiplicity.fArray < (2.536e+7 + TMath::Sqrt(1.574e14+1.390e12*(lsv_cluster_fixed_width_charge.fArray-14.40)**2))";
      odtree->Draw("lsv_cluster_fixed_width_charge.fArray>>FullSpectrum", TimingCut && MultiCut);
      //  FullSpectrum->Rebin(4);
      //  FullSpectrum->Sumw2();
      FullSpectrum->Scale(1.0/nEntries);
      */
      //    }

  if(tag==1){
    // string MCDataFile = split->GetMCinputdir() + "PulseSplitMCEnergy.root";
    // TFile *MCFile = new TFile(MCDataFile.c_str());
    split->Recon_DataAnalysis();
    string RealDataFile = split->GetRealinputdir() +"PulseSplitRealEnergyUAr_mtree_May2.root";
    //  string RealDataFile = split->GetRealinputdir() +"PulseSplitRealEnergy_MoreBins.root";
    cout<<"Fitting "<<RealDataFile<<endl;
    TFile *RealFile = new TFile(RealDataFile.c_str());    
    TH1F *FullSpectrum = (TH1F*) RealFile->Get("C14Spectrum");
    FullSpectrum->SetTitle("");
    TLegend *leg = new TLegend(0.4,0.75,0.7,0.9);
    TCanvas *c1 = new TCanvas("c1", "Full Spectrum with C14 Fit", 1000, 400);
    c1->SetLogy();
    c1->cd();
    FullSpectrum->Draw();
    double startrange=25;
    double endrange = 150;
    int startbin = FullSpectrum->FindBin(startrange);
    int endbin = FullSpectrum->FindBin(endrange);
    Double_t integralsum = FullSpectrum->Integral(startbin,endbin);
    const Int_t parnums =10;
    /*    const Int_t isotopes = 6;
	  string Source[] = {"C14","C14_1st","C14_2nd","Co60","Co57","K40","Th232","U238"};
    Double_t isotope_activity[isotopes] = {223268845,223268845,223268845,99473576,24773624,47584434};
    Double_t isotope_fraction[isotopes];
    split->Isotope_Fraction(isotopes,isotope_activity,isotope_fraction);
    Double_t integral_sum[isotopes];
    for(int i=0; i<isotopes; i++)
      integral_sum[i] = isotope_fraction[i]*integralsum;
      Double_t kB               = 0.00835;
    */
    Double_t new_par[parnums];    
    Double_t binw             = FullSpectrum->GetBinWidth(1);
    /*
    Double_t ly_mean          = 0.55;
    Double_t spe_var          = 0.14;
    Double_t ly_var           = 0.00642;
    Double_t constant         = 0.;
    Double_t basel_mean       = 0.386;
    Double_t basel_var        = 1.;
    Double_t threshold        = 1;
    Double_t window           = 5.78e-7; // 2e-7;
    */
    Double_t ly_mean          = 0.55;
    Double_t spe_var          = 0.14;
    Double_t ly_var           = 0.00642;
    Double_t constant         = 1;
    Double_t basel_mean       = 0.386;
    Double_t basel_var        = 100.;
    Double_t threshold        = 10.;
    Double_t window           = 5.78e-7; // 2e-7;
    
    Double_t parvaluesamples[]={binw,ly_mean,spe_var,ly_var,constant,basel_mean,basel_var,threshold,window,integralsum};
				//integralsum,integralsum,integralsum,integralsum,integralsum,integralsum,integralsum};
				//			integral_sum[0],integral_sum[1],integral_sum[2],integral_sum[3],integral_sum[4],integral_sum[5]};
    string parnamesamples[]={"Bin Width","LY Mean [PE/keV]","Rel SPE Var","Rel LY Var",
			     "Constant","Baseline Mean [p.e]","Baseline Var","Threshold","Window Width [s]","C14 Rate [Hz]"};
			     //		     "C14 Amplitude","Co60 Amplitude","Co57 Amplitude","K40 Amplitude","Th232 Amplitude","U238 Amplitude"};
    Double_t parlowlimits[parnums] = {0,0,0,0,0,0,0,0,0,0};
    Double_t paruplimits[]  = {10,1.,0.14,100,500,500,500,500,1e-6,integralsum*5};
    
    TF1* Fit_C14 = new TF1("Fit_C14",split,&c14pulse::C14Fit, startrange, endrange, parnums,"c14pulse","C14Fit");        
    TF1* Fit_C14_1st = new TF1("Fit_C14_1st",split,&c14pulse::C14Fit_1st , startrange, endrange, parnums,"c14pulse","C14Fit_1st"); 
    TF1* Fit_C14_2nd = new TF1("Fit_C14_2nd", split,&c14pulse::C14Fit_2nd, startrange, endrange, parnums,"c14pulse","C14Fit_2nd"); 
    TF1* Fit_Total  = new TF1("Fit_Total",split,&c14pulse::TotalFit,startrange,endrange,parnums,"c14pulse","TotalFit");
    //   FitC14 = new TF1("FitC14",C14Fit,start[0],end[0],parnums); //total C14 fit function
    /*    TF1* Fit_Co60  = new TF1("Fit_Co60",split,&c14pulse::Co60Fit,startrange,endrange,parnums,"c14pulse","Co60Fit");
    TF1* Fit_Co57  = new TF1("Fit_Co57",split,&c14pulse::Co57Fit,startrange,endrange,parnums,"c14pulse","Co57Fit");
    TF1* Fit_K40   = new TF1("Fit_K40",split,&c14pulse::K40Fit,startrange,endrange,parnums,"c14pulse","K40Fit");
    //  TF1* Fit_Tl208 = new TF1("Fit_Tl208",split,&c14pulse::Tl208Fit,startrange,endrange,parnums,"c14pulse","Tl208Fit");
    TF1* Fit_Th232 = new TF1("Fit_Th232",split,&c14pulse::Th232Fit,startrange,endrange,parnums,"c14pulse","Th232Fit");
    TF1* Fit_U238  = new TF1("Fit_U238",split,&c14pulse::U238Fit,startrange,endrange,parnums,"c14pulse","U238Fit");
    */
    TObjArray Fitlist(0);
    Fitlist.Add(Fit_C14);
    Fitlist.Add(Fit_C14_1st);
    Fitlist.Add(Fit_C14_2nd);
    //    Fitlist.Add(Fit_Co60);
    //    Fitlist.Add(Fit_Co57);
    //    Fitlist.Add(Fit_K40);
    // Fitlist.Add(Fit_Tl208);
    //    Fitlist.Add(Fit_Th232);
    //    Fitlist.Add(Fit_U238);
    //  Fitlist.Add(Fit_Total);
    
    vector<int> linecolor = split->Colors();
    for(int i=0; i<parnums; i++)
      {
	Fit_Total->SetParName(i,parnamesamples[i].c_str());
	if(i==3 || i==5 || i==8 || i == 0 || i==2) Fit_Total->FixParameter(i,parvaluesamples[i]);
	else{
	  Fit_Total->SetParameter(i, parvaluesamples[i]);
	  Fit_Total->SetParLimits(i,parlowlimits[i],paruplimits[i]);
	}
      }
    Fit_Total->SetNpx(200);
    //  FitTotal->SetLineColor(kRed);
    Fit_Total->SetLineColor(linecolor.at(Fitlist.GetEntries()+1));  

    FullSpectrum->Fit(Fit_Total,"RVEM");    
    Fit_Total->GetParameters(new_par);
    //    Double_t* new_par_error; 
    //    new_par_error = Fit_Total->GetParErrors();

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

    leg->AddEntry(Fit_Total,"Total Fit","l");
    leg->Draw();

    Double_t chi2 = Fit_Total->GetChisquare();
    Int_t ndf = Fit_Total->GetNDF();
    cout<<"chi2/ndf= "<<chi2/(1.0*ndf)<<endl;
    string outdir=split->Getoutputdir();
    string output = outdir + "C14Fit_"+Time+".root";
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




