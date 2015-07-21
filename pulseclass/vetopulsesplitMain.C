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
#include "vetopulsesplit.hh"
#include "./odselector/odselector.h"
#include "./DSTtreeSelector/DSTtreeSelector.h"

using namespace std;
//#define REBINS

TRint* theApp;
string Time = "Jun10";
string Tag = "v4t";
string dirname="/darkside/users/hqian/pulse_splitter/";
#define fieldon

void Process(TChain* chain, TString label){
  TString option = label;
  chain->SetProof();
  TProof* pr = TProof::Open("workers=2");
  pr->SetParameter("PROOF_Packetizer","TPacketizer");
  pr->SetParameter("PROOF_MaxSlavesPerNode",8);
  chain->Process("./odselector/odselector.C+",option.Data(),-1,0);
}

void DST_Process(TChain* chain, TString label){
  TString option = label;
  chain->SetProof();
  TProof* pr = TProof::Open("workers=2");
  pr->SetParameter("PROOF_Packetizer","TPacketizer");
  pr->SetParameter("PROOF_MaxSlavesPerNode",8);
  chain->Process("./DSTtreeSelector/DSTtreeSelector.C+",option.Data(),-1,0);
}

#ifndef __CINT__
int main(int argc, char** argv){
  theApp = new TRint("App",&argc,argv,NULL,0);
  int start, end, tag;
  bool start_fit;
  if(theApp->Argc() == 1)
    tag = 1;
  else if(theApp->Argc() == 2)
    {
      tag = 2;
      start_fit = atoi(theApp->Argv(1))>0 ? true : false;
    }
  else if(theApp->Argc() == 3)
    {
      tag=3;
      start = atoi(theApp->Argv(1));
      end = atoi(theApp->Argv(2));
    }
  else if(theApp->Argc() == 4)
    {
      tag= atoi(theApp->Argv(3));
      start = atoi(theApp->Argv(1));
      end = atoi(theApp->Argv(2));
    }
  else{
    cout<<"Usage: ./reconMain startfile endfile "<<endl;
    cout<<"Usage: ./reconMain "<<endl;
    return 0;
  }

  vetopulsesplit *split = new vetopulsesplit;
  //  split->SetMCinputdir(dirname+"SDMCData/");
  split->SetMCinputdir(dirname+"50daysMCData/");
  split->SetRealinputdir(dirname+"VetoData/");
  split->Setoutputdir(dirname);  

  //  string realoutfile = "PulseSplitRealEnergy_MoreBins_noerr.root";
#ifdef fieldon
  string realoutfile = "PulseSplitRealEnergyUAr_mtree_May5.root";
#else
  string realoutfile = "PulseSplitRealEnergyUAr_nullfield_mtree_May5.root";
#endif

  //  string mcfile = "PulseSplitMCEnergy_"+Tag+"_May22.root";
  string mcfile = "PulseSplitMCEnergy_"+Tag+"_"+Time+".root";
  cout<<mcfile<<endl; 
  split->SetMCFile(mcfile);  
  split->SetIsotopes();

  if(tag==3)
    {
      //      split->ODTree_DataAnalysis(start,end);      
      TChain *odtree =new TChain("odtree");
      split->Readdatafile(odtree,start,end);
      Int_t nEntries = odtree->GetEntries();
      cout<<"nEntries= "<<nEntries<<endl;
      string label = split->GetRealinputdir() + realoutfile;      
      Process(odtree,label);
    }
  else if(tag==4){
      TChain *DSTtree =new TChain("DSTtree");
      split->DST_Readdatafile(DSTtree,start,end);
      Int_t nEntries = DSTtree->GetEntries();
      cout<<"nEntries= "<<nEntries<<endl;
      // string label = split->GetRealinputdir() + realoutfile;      
      string label = realoutfile;      
      DST_Process(DSTtree,label);
    }
  else if(tag==1){
    split->Recon_DataAnalysis();
  }
  else if(tag==2){
    //    split->Recon_DataAnalysis();
    string MCDataFile = split->GetMCinputdir() + mcfile;
    cout<<MCDataFile<<endl;
    TFile *MCFile = new TFile(MCDataFile.c_str());
    vector<TH1F*> mcplot;
    TObjArray MClist(0);
    vector<string> Source = split->GetIsotopes();
    for(size_t i=0; i<Source.size(); i++)
      {
	mcplot.push_back((TH1F*)MCFile->Get(Form("%s_EnergyMC",Source.at(i).c_str())));
	MClist.Add(mcplot.at(i));
      }
    split->SetEnergyMC(mcplot);
    TNtuple *MCValue_ntuple = (TNtuple*) MCFile->Get("MCValue_ntuple");
    vector<double> isotope_activity;
    float fraction;
    MCValue_ntuple->SetBranchAddress("fraction",&fraction);
    for(int i=0; i<MCValue_ntuple->GetEntries(); i++)
      {
	MCValue_ntuple->GetEntry(i);
	isotope_activity.push_back(fraction);	
      }
    string RealDataFile = split->GetRealinputdir() + realoutfile;
    cout<<"Fitting "<<RealDataFile<<endl;
    TFile *RealFile = new TFile(RealDataFile.c_str());    
    TH1F *FullSpectrum = (TH1F*) RealFile->Get("FullSpectrum");
    int rebins = 5;
    FullSpectrum->Rebin(rebins);    
    FullSpectrum->Sumw2();  
#ifdef fieldon
    double livetime = 449177; //6103.64; //[s]
#else 
    double livetime = 253431; //[s]
#endif
    FullSpectrum->Scale(1./livetime);      
    FullSpectrum->SetTitle("");

    TLegend *leg = new TLegend(0.5,0.7,0.8,1.0);
    TCanvas *c1 = new TCanvas("c1", "Full Spectrum with Fit", 1000, 400);
    c1->SetLogy();
    c1->cd();
    FullSpectrum->Draw();
    double startrange= 0;
    double endrange = 2000;
    int Bins = FullSpectrum->GetNbinsX();
    int startbin = FullSpectrum->FindBin(startrange);
    int endbin = FullSpectrum->FindBin(endrange);
    Double_t integralsum = FullSpectrum->Integral(startbin,endbin);
    cout<<"total integralsum= "<<integralsum<<endl;
    
    //    vector<string> Source = split->GetIsotopes();
    int isotopes = static_cast<int> (Source.size());
    const Int_t parnums =  isotopes + 9;

    //    vector<double> isotope_activity = split->GetFraction();
    Double_t isotope_fraction[isotopes];
    split->Isotope_Fraction(isotopes,isotope_activity,isotope_fraction);
    Double_t integral_sum[isotopes];
    for(int i=0; i<isotopes; i++)
      {
	integral_sum[i] = isotope_fraction[i]*integralsum;    
	//	integral_sum[i] = integralsum;    
	cout<<"Source[i]= "<<Source.at(i)<<"\t isotope_fraction[i] "<<isotope_fraction[i]<<
	  "\t integral_sum="<<integral_sum[i]<<endl;
      }
    Double_t new_par[parnums];
    //  Double_t kB               = 0.00835;     
    Double_t binw             = FullSpectrum->GetBinWidth(1);
    Double_t ly_mean          = 0.6058;  //0.63;
    Double_t spe_var          = 0.14;
    Double_t ly_var           = 0.0043;  // 0.00642;
    Double_t constant         = 0.;
    Double_t basel_mean       = 0.386;
    Double_t basel_var        = 100; //1.;
    Double_t threshold        = 0; //1;
    Double_t window           = 5.78e-7;   // 2e-7;
       
    Double_t parvaluesamples[]={binw,ly_mean,spe_var,ly_var,constant,basel_mean,basel_var,threshold,window,
				//integralsum,integralsum,integralsum,integralsum,integralsum,integralsum,integralsum};
				integral_sum[0],integral_sum[1],integral_sum[2],integral_sum[3],
				integral_sum[4],integral_sum[5],0,integral_sum[7],integral_sum[8]*10,integral_sum[9]*10};
    string parnamesamples[]={"Bin Width","LY Mean [PE/keV]","Rel SPE Var","Rel LY Var",
			     "Constant","Baseline Mean [p.e]","Baseline Var","Threshold","Window Width [s]",
			     "C14 Rate[Bq]","Co60 Rate[Bq]","Co57 Rate[Bq]","K40 Rate[Bq]",
			     "Th232Upper Rate[Bq]","Th232Lower Rate[Bq]","U235Upper Rate[Bq]","U235Lower Rate[Bq]","U238Upper Rate[Bq]","U238Lower Rate[Bq]"};
    // Double_t parlowlimits[parnums] = {0,0,0,0,-100,-10,-100,0,0,0,0,0,0,0,0};
    // Double_t paruplimits[]  = {10,1.,0.14,100,100,100,100,100,1e-6,integralsum*1,integralsum*2,integralsum*2,integralsum*2,integralsum*2,integralsum*2}; 

    Double_t parlowlimits[] = {0,0,0,0,0,0,0,0,0,
			       // 1,1,1,1,1,1,1};
			       0,0,0,0,0,0,0,0,0,0};
    Double_t paruplimits[]  = {20,1,1,1,10,50,100,100,1e-6,
			       integralsum*20,integralsum*1,integralsum*20,integralsum*20,
			       integralsum*100,integralsum*100,integralsum*100/21.5,integralsum*100,integralsum*100,integralsum*100};
    
    TF1* Fit_C14 = new TF1("Fit_C14",split,&vetopulsesplit::C14Fit, startrange, endrange, parnums,"vetopulsesplit","C14Fit");        
    //    TF1* Fit_C14_1st = new TF1("Fit_C14_1st",split,&vetopulsesplit::C14Fit_1st , startrange, endrange, parnums,"vetopulsesplit","C14Fit_1st"); 
    // TF1* Fit_C14_2nd = new TF1("Fit_C14_2nd", split,&vetopulsesplit::C14Fit_2nd, startrange, endrange, parnums,"vetopulsesplit","C14Fit_2nd"); 
    //   FitC14 = new TF1("FitC14",C14Fit,start[0],end[0],parnums); //total C14 fit function
    TF1* Fit_Co60  = new TF1("Fit_Co60",split,&vetopulsesplit::Co60Fit,startrange,endrange,parnums,"vetopulsesplit","Co60Fit");
    TF1* Fit_Co57  = new TF1("Fit_Co57",split,&vetopulsesplit::Co57Fit,startrange,endrange,parnums,"vetopulsesplit","Co57Fit");
    TF1* Fit_K40   = new TF1("Fit_K40",split,&vetopulsesplit::K40Fit,startrange,endrange,parnums,"vetopulsesplit","K40Fit");
    //  TF1* Fit_Tl208 = new TF1("Fit_Tl208",split,&vetopulsesplit::Tl208Fit,startrange,endrange,parnums,"vetopulsesplit","Tl208Fit");
    TF1* Fit_Th232 = new TF1("Fit_Th232",split,&vetopulsesplit::Th232Fit,startrange,endrange,parnums,"vetopulsesplit","Th232Fit");
    TF1* Fit_Th232Lower = new TF1("Fit_Th232Lower",split,&vetopulsesplit::Th232LowerFit,startrange,endrange,parnums,"vetopulsesplit","Th232LowerFit");
    TF1* Fit_U235  = new TF1("Fit_U235",split,&vetopulsesplit::U235Fit,startrange,endrange,parnums,"vetopulsesplit","U235Fit");
    TF1* Fit_U235Lower  = new TF1("Fit_U235Lower",split,&vetopulsesplit::U235LowerFit,startrange,endrange,parnums,"vetopulsesplit","U235LowerFit");
    TF1* Fit_U238Upper  = new TF1("Fit_U238Upper",split,&vetopulsesplit::U238UpperFit,startrange,endrange,parnums,"vetopulsesplit","U238UpperFit");
    TF1* Fit_U238Lower  = new TF1("Fit_U238Lower",split,&vetopulsesplit::U238LowerFit,startrange,endrange,parnums,"vetopulsesplit","U238LowerFit");
    //    TF1* Fit_U238  = new TF1("Fit_U238",split,&vetopulsesplit::U238Fit,startrange,endrange,parnums,"vetopulsesplit","U238Fit");
    //TF1* Fit_Rn222  = new TF1("Fit_Rn222",split,&vetopulsesplit::Rn222Fit,startrange,endrange,parnums,"vetopulsesplit","Rn222Fit");
    //    TF1* Fit_K42  = new TF1("Fit_K42",split,&vetopulsesplit::K42Fit,startrange,endrange,parnums,"vetopulsesplit","K42Fit");
    //    TF1* Fit_Rn222  = new TF1("Fit_Rn222",split,&vetopulsesplit::Rn222Fit,startrange,endrange,parnums,"vetopulsesplit","Rn222Fit");
    TF1* Fit_Total  = new TF1("Fit_Total",split,&vetopulsesplit::TotalFit,startrange,endrange,parnums,"vetopulsesplit","TotalFit");

    TObjArray Fitlist(0);
    Fitlist.Add(Fit_C14);
    //    Fitlist.Add(Fit_C14_1st);
    //    Fitlist.Add(Fit_C14_2nd);
    Fitlist.Add(Fit_Co60);
    Fitlist.Add(Fit_Co57);
    Fitlist.Add(Fit_K40);
    // Fitlist.Add(Fit_Tl208);
    Fitlist.Add(Fit_Th232);
    Fitlist.Add(Fit_Th232Lower);
    Fitlist.Add(Fit_U235);
    Fitlist.Add(Fit_U235Lower);
    Fitlist.Add(Fit_U238Upper);
    Fitlist.Add(Fit_U238Lower);
    //    Fitlist.Add(Fit_U238);
    //Fitlist.Add(Fit_Rn222);
    //  Fitlist.Add(Fit_K42);

    vector<int> linecolor = split->Colors();
    for(int i=0; i<parnums; i++)
      {
       	Fit_Total->SetParName(i,parnamesamples[i].c_str());
       	if(i==0||i==2||i==3||i==4||i==7||i==8||i==5||i==15)	
	  Fit_Total->FixParameter(i,parvaluesamples[i]);
	else 
	  { 
	    Fit_Total->SetParameter(i, parvaluesamples[i]);
	    Fit_Total->SetParLimits(i,parlowlimits[i],paruplimits[i]);
	  }
      }
    Fit_Total->SetNpx(2000);
    Fit_Total->SetLineColor(linecolor.at(12));  

    if(start_fit)
      {
	//FullSpectrum->Fit(Fit_Total,"RVWLM");
        FullSpectrum->Fit(Fit_Total,"RVEM");
      }
    Fit_Total->GetParameters(new_par);
    for(Int_t i=0; i<Fitlist.GetEntries(); i++)
      {
	TF1* tempfit = dynamic_cast<TF1*>(Fitlist.At(i));
	tempfit->SetLineColor(linecolor.at(i));
	tempfit->SetParameters(new_par);
	tempfit->Draw("SAME");
	leg->AddEntry(tempfit,tempfit->GetName(),"l");
	cout<<i<<"\t"<<tempfit->GetName()<<"  "<<Fitlist.GetSize()<<" "<<Fitlist.GetEntries()<<endl;
      }
    Fit_Total->Draw("same");
    leg->AddEntry(Fit_Total,"Total Fit","l");
    leg->Draw();

    Double_t chi2 = Fit_Total->GetChisquare();
    Int_t ndf = Fit_Total->GetNDF();
    cout<<"chi2/ndf= "<<chi2/(1.0*ndf)<<"\t Probability="<<Fit_Total->GetProb()<<endl;
    string outdir=split->Getoutputdir();
#ifdef fieldon
    string output = outdir + "Pulse_Splitter_"+Tag+"_fieldon_"+Time+".root";
#else
    string output = outdir + "Pulse_Splitter_"+Tag+"_fieldoff_"+Time+".root";
#endif
    TFile outfile(output.c_str(), "RECREATE");
    //    split->GetHistArray().Write();
    MClist.Write();
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










    /*
    //#ifdef REBINS
    TNtuple *FullSpectrum_ntuple = (TNtuple*) RealFile->Get("FullSpectrum_ntuple");
    const int seperated = 2;
    int seperated_range[seperated] = {1000,2000}; 
    int seperated_binw[seperated] = {5,10};
    int seperated_bin[seperated];    
    int nbinsx_temp=0;
    for(int i=0; i<seperated; i++)
      {
	if(i==0)
	  nbinsx_temp += seperated_range[i]/seperated_binw[i];
	else
	  nbinsx_temp += (seperated_range[i]-seperated_range[i-1])/seperated_binw[i];
	seperated_bin[i] = nbinsx_temp;
      }
    const int nbinsx = nbinsx_temp; 
    double xbins[nbinsx+1];     
    int k=0;
    for(int i=0; i<seperated; i++)
      {
	while(k<seperated_bin[i])	      
	  {
	    if(i==0)
	      xbins[k] = seperated_binw[i]*k;
	    else  xbins[k] = seperated_binw[i]*(k-seperated_bin[i-1])+seperated_range[i-1]; 
	    cout<<"k= "<<k<<"\t xbins[k]="<<xbins[k]<<endl;
	    k++;
	  }
      }	
    xbins[nbinsx]=seperated_range[seperated-1];
    TH1F *FullSpectrum = new TH1F("FullSpectrum_new","",nbinsx,xbins);   
    Float_t         charge;
    FullSpectrum_ntuple->SetBranchAddress("charge", &charge);    
    for(int j=0; j<FullSpectrum_ntuple->GetEntries(); j++)
      {
	FullSpectrum_ntuple->GetEntry(j);
	for(int i=0; i<seperated; i++)
	  {
	    if(charge<=seperated_range[i])
	      FullSpectrum->Fill(charge,1./seperated_binw[i]);
	    else continue;
	  }
      }
    FullSpectrum->Sumw2();
    */
