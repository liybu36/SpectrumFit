#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooFormulaVar.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "RooFitResult.h"

#include "TApplication.h"
#include "TRint.h"
#include "TProof.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TChain.h"
#include "TFile.h"
#include "TNtuple.h"
#include "../vetopulsesplit.hh"
#include "../odselector/odselector.h"
#include "../DSTtreeSelector/DSTtreeSelector.h"

using namespace std;
using namespace RooFit;

TRint* theApp;
string Time = "Apr10PM";
string Tag = "v4";
//string dirname="/darkside/users/hqian/pulse_splitter/";
string dirname="/Users/hqian/Results/Apr27/";

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
  split->SetMCinputdir(dirname+"MCData/");
  split->SetRealinputdir(dirname+"VetoData/");
  split->Setoutputdir(dirname);

  //  string realoutfile = "PulseSplitRealEnergy_MoreBins_noerr.root";
  //  string realoutfile = "PulseSplitRealEnergyUAr_nullfield.root";
  //  string realoutfile = "PulseSplitRealEnergyUAr.root";
  string realoutfile = "PulseSplitRealEnergyUAr_tree.root";

  string mcfile = "PulseSplitMCEnergy_"+Tag+".root";
  cout<<mcfile<<endl;
  split->SetMCFile(mcfile);
  split->SetIsotopes();
  
  if(tag==1){
    split->Recon_DataAnalysis();
  }
  else if(tag==3)
    {
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
    string label = split->GetRealinputdir() + realoutfile;
    DST_Process(DSTtree,label);
  }
  else if(tag==2){
    string MCDataFile = split->GetMCinputdir() + mcfile;
    cout<<MCDataFile<<endl;
    TFile *MCFile = new TFile(MCDataFile.c_str());
    vector<TH1F*> mcplot;
    TObjArray MClist(0);
    vector<string> Source = split->GetIsotopes();
    vector<int> NBins = split->GetNBins();
    vector<TString> hname;
    vector<TString> htitle;    
    for(size_t i=0; i<Source.size(); i++)
      {
	hname.push_back(Form("%s_MCEnergy",Source.at(i).c_str()));
	htitle.push_back(Form("%s MC Energy",Source.at(i).c_str()));
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

    RooRealVar charge("charge","charge",0,2000);

    TTree* FullSpectrum_tree = (TTree*) RealFile->Get("lsvtree");
    RooDataSet data("data","Real Data",FullSpectrum_tree,charge);
    data.Print();
    
    /*
    TH1F *FullSpectrum = (TH1F*) RealFile->Get("FullSpectrum");
    RooDataHist data("data","Real Data",RooArgSet(charge),Import(*FullSpectrum));
    */
    RooRealVar binw("binw","Bin Width",1,0,10);
    RooRealVar ly_mean("ly_mean","LY Mean[PE/keV]",0.63,0,1);
    RooRealVar spe_var("spe_var","Rel SPE Var",0.14,0,1);
    RooRealVar ly_var("ly_var","Rel LY Var",0.00642,0,1);
    RooRealVar constant("constant","Constant",0,0,10);
    RooRealVar basel_mean("basel_mean","Baseline Mean[PE]",0.386,0,500);
    RooRealVar basel_var("basel_var","Baseline Var",100,0,500);
    RooRealVar threshold("threshold","Threshold",10,0,500);
    RooRealVar window("window","Window Width[s]",5.78e-7,0,1.e-6);

    //    binw.SetConstant(kTRUE);
    //    constant.SetConstant(kTRUE);
    //    window.SetConstant(kTRUE);
    
    vector<RooRealVar*> mc_x;
    vector<RooDataHist*> mc_hist;
    vector<RooHistPdf*> mc_pdf;
    vector<RooFormulaVar*> mc_gaus_mean;
    vector<RooFormulaVar*> mc_gaus_var;
    vector<RooGaussian*> mc_gauss;
    vector<RooFFTConvPdf*> mc_convpdf;
    vector<RooRealVar*> mc_fraction;    
    RooArgList pdflist("pdflist");
    RooArgList fraclist("fraclist");
    
    for(size_t i=0; i<Source.size(); i++)
      {
	TString temp_xname,temp_xtitle,temp_pdfname,temp_pdftitle,temp_gaus_mean_name,temp_gaus_var_name,
	  temp_gauss_name,temp_gauss_title,temp_convpdfname,temp_convpdftitle;
	temp_xname.Form("%s_x",Source.at(i).c_str());
	temp_xtitle.Form("%s x",Source.at(i).c_str());
	RooRealVar *temp_x = new RooRealVar(temp_xname,temp_xtitle,0,NBins.at(i));
	//RooRealVar temp_x(temp_xname,temp_xtitle,0,NBins.at(i));
	mc_x.push_back(temp_x);
	RooDataHist *temp_hist = new RooDataHist(hname.at(i),htitle.at(i),RooArgSet(*(mc_x.at(i))),mcplot.at(i));    
	mc_hist.push_back(temp_hist);
	temp_pdfname.Form("%s_pdf",Source.at(i).c_str());
	temp_pdftitle.Form("%s pdf",Source.at(i).c_str());	
	RooHistPdf *temp_pdf = new RooHistPdf(temp_pdfname,temp_pdftitle,*(mc_x.at(i)),*(mc_hist.at(i)),0);
   	mc_pdf.push_back(temp_pdf);
	
	temp_gaus_mean_name.Form("%s_gaus_mean",Source.at(i).c_str());
        temp_gaus_var_name.Form("%s_gaus_var",Source.at(i).c_str());
	RooFormulaVar *temp_gaus_mean = new RooFormulaVar(temp_gaus_mean_name,"@0*@1+@2",RooArgList(*(mc_x.at(i)),ly_mean,basel_mean));
	RooFormulaVar *temp_gaus_var = new RooFormulaVar(temp_gaus_var_name,"@0+(1+@1)*@2*@3+@4*TMath::Power(@2*@3,2)", 
							 RooArgList(basel_var,spe_var,*(mc_x.at(i)),ly_mean,ly_var));    
		
	//"basel_var+(1+spe_var)*energy*ly_mean+ly_var*TMath::Power(energy*ly_mean,2)",
	//	RooFormulaVar temp_gaus_mean(temp_gaus_mean_name,"@0*@1+@2",RooArgList(mc_x.at(i),ly_mean,basel_mean));
	//	RooFormulaVar temp_gaus_var(temp_gaus_var_name,"@0+(1+@1)*@2*@3+@4*TMath::Power(@2*@3,2)", 			      
	//				    RooArgList(basel_var,spe_var,mc_x.at(i),ly_mean,ly_var));    

	mc_gaus_mean.push_back(temp_gaus_mean);
	mc_gaus_var.push_back(temp_gaus_var);
        temp_gauss_name.Form("%s_gauss",Source.at(i).c_str());
	temp_gauss_title.Form("%s Gaussian",Source.at(i).c_str());
	RooGaussian *temp_gauss = new RooGaussian(temp_gauss_name,temp_gauss_title,charge,*(mc_gaus_mean.at(i)),*(mc_gaus_var.at(i)));
	mc_gauss.push_back(temp_gauss);

	temp_convpdfname.Form("%s_convpdf",Source.at(i).c_str());
	temp_convpdftitle.Form("%s Signal & Response convolution",Source.at(i).c_str());	
	RooFFTConvPdf *temp_convpdf = new RooFFTConvPdf(temp_convpdfname,temp_convpdftitle,charge,*(mc_pdf.at(i)),*(mc_gauss.at(i)));
	mc_convpdf.push_back(temp_convpdf);
	pdflist.add(*(mc_convpdf.at(i)));
	
	RooRealVar *temp_fraction = new RooRealVar(Form("%s_frac",Source.at(i).c_str()),Form("%s fraction",Source.at(i).c_str()),isotope_activity.at(i),0,1);
        // RooRealVar temp_fraction(Form("%s_frac",Source.at(i).c_str()),Form("%s fraction",Source.at(i).c_str()),isotope_activity.at(i),0,1);
	mc_fraction.push_back(temp_fraction);
	fraclist.add(*(mc_fraction.at(i)));

      }
    //    Bool_t IsRecursive = true;
    // RooAddPdf model("model","Total Fit",pdflist,fraclist,IsRecursive);	
    RooAddPdf model("model","Total Fit",pdflist,fraclist);	

    if(start_fit){
      RooFitResult *r = model.fitTo(data,SumW2Error(kTRUE),Save());    
      r->Print();
    }
    
    RooPlot* frame = charge.frame();
    data.plotOn(frame);
    model.plotOn(frame);
    
    TCanvas *c1 = new TCanvas("c1","Full Spectrum Fit",1000,600);
    frame->Draw();

  }
  return 1;   
  
}
#endif /*__CINT__*/











    /*    
    RooRealVar c14_rate("c14_rate","C14 Rate[Hz]",sum.at(0),0,integralsum);
    RooRealVar co60_rate("co60_rate","Co60 Rate[Hz]",sum.at(1),0,integralsum);
    RooRealVar co57_rate("co57_rate","Co57 Rate[Hz]",sum.at(2),0,integralsum);
    RooRealVar k40_rate("k40_rate","K40 Rate[Hz]",sum.at(3),0,integralsum);
    RooRealVar th232_rate("th232_rate","Th232 Rate[Hz]",sum.at(4),0,integralsum);
    RooRealVar u238_rate("u238_rate","U238Upper Rate[Hz]",sum.at(5),0,integralsum);
    RooRealVar u238lower_rate("u238lower_rate","U238Lower Rate[Hz]",sum.at(6),0,integralsum);
    RooRealVar u235_rate("u235_rate","U235 Rate[Hz]",sum.at(7),0,integralsum);
        
    RooRealVar x("x","x",0,2000);
    RooFormulaVar gaus_mean("gaus_mean","x*ly_mean+basel_mean",RooArgList(x,ly_mean,basel_mean));
    RooFormulaVar gaus_var("gaus_var","@0+(1+@1)*@2*@3+@4*TMath::Power(@2*@3,2)", 
			   //"basel_var+(1+spe_var)*energy*ly_mean+ly_var*TMath::Power(energy*ly_mean,2)",
			   RooArgList(basel_var,spe_var,x,ly_mean,ly_var));    
    RooGaussian gauss("gauss","gauss",x,gaus_mean,gaus_var);

    x.SetBins("fft",2000*2);

    RooRealVar c14_x("c14_x","c14_x",0,NBins[0]);
    RooDataHist c14_hist(hname.at(0),htitle.at(0),RooArgSet(c14_x),mcplot.at(0));    
    RooHistPdf c14_pdf("c14_pdf","c14_pdf",c14_x,c14_hist,0) ;
    RooFFTConvPdf conv("conv","Siganl (x) Response",x,c14_pdf,gauss);

    RooRealVar co60_x("co60_x","co60_x",0,NBins[1]);
    RooDataHist co60_hist(hname.at(1),htitle.at(1),RooArgSet(co60_x),mcplot.at(1));
    RooHistPdf co60_pdf("co60_pdf","co60_pdf",co60_x,co60_hist,0);
    
    RooRealVar eff("eff","signal fraction",0,1);
    RooAddPdf model("model","Total Fit",RooArgList(c14_pdf,co60_pdf),eff);
    */
