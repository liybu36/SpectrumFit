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
#include "RooLandau.h"

//#include "TApplication.h"
//#include "TRint.h"
#include "TProof.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TChain.h"
#include "TFile.h"
#include "TNtuple.h"

using namespace std;
using namespace RooFit;

//TRint* theApp;
string Time = "Apr10PM";
string Tag = "v4";

vector<string> GetIsotopes()
{
  vector<string> Source;
  
  Source.push_back("C14");
  Source.push_back("Co60");
  Source.push_back("Co57");
  Source.push_back("K40");
  Source.push_back("Th232");
  Source.push_back("U238");
  Source.push_back("U238Lower");
  Source.push_back("U235");

  return Source;
}
 
vector<int> GetNBins()
{
  vector<int> NBins;

  NBins.push_back(160);
  NBins.push_back(2600);
  NBins.push_back(800);
  NBins.push_back(1600);
  NBins.push_back(3500);
  NBins.push_back(3500);
  NBins.push_back(3500);
  NBins.push_back(1000);

  return NBins;
}

vector<int> Colors()
{
  vector<int> color;
  color.push_back(TColor::GetColor("#FF0000")); //red         0
  color.push_back(TColor::GetColor("#5A1DE8")); //violet      1
  color.push_back(TColor::GetColor("#000000"));
  //  color.push_back(TColor::GetColor("#E8A60C")); //brown       2
  color.push_back(TColor::GetColor("#F73CFF")); //pink        5
  color.push_back(TColor::GetColor("#1CFFDF")); //low green   7
  color.push_back(TColor::GetColor("#1485CC")); //blue        4
  color.push_back(TColor::GetColor("#FF791F")); //orange      6
  color.push_back(TColor::GetColor("#AF2FCC")); //dark pink   8
  color.push_back(TColor::GetColor("#E8A60C"));
  color.push_back(TColor::GetColor("#59FF49")); //green       9

  return color;
}

int vetospectrumroofit()
{
  bool start_fit = true;
  //#define single_fit
  size_t iter = 0;

  string dirname="/darkside/users/hqian/pulse_splitter/";
  //string dirname="/Users/hqian/Results/Apr27/";
  string mcdir = dirname+"MCData/";
  string realdir = dirname+"VetoData/";
    
  //  string realoutfile = "PulseSplitRealEnergy_MoreBins_noerr.root";
  //  string realoutfile = "PulseSplitRealEnergyUAr_nullfield.root";
  //  string realoutfile = "PulseSplitRealEnergyUAr.root";
  string realoutfile = "PulseSplitRealEnergyUAr_tree.root";
  string mcfile = "PulseSplitMCEnergy_"+Tag+".root";

  string MCDataFile = mcdir + mcfile;
  cout<<MCDataFile<<endl;
  TFile *MCFile = new TFile(MCDataFile.c_str());
  vector<TH1F*> mcplot;
  TObjArray MClist(0);
  vector<string> Source = GetIsotopes();
  vector<int> NBins = GetNBins();
  vector<int> linecolor = Colors();
  vector<TString> hname;
  vector<TString> htitle;    
  for(size_t i=0; i<Source.size(); i++)
    {
      hname.push_back(Form("%s_MCEnergy",Source.at(i).c_str()));
      htitle.push_back(Form("%s MC Energy",Source.at(i).c_str()));
      mcplot.push_back((TH1F*)MCFile->Get(Form("%s_EnergyMC",Source.at(i).c_str())));
      MClist.Add(mcplot.at(i));	
    }

  TNtuple *MCValue_ntuple = (TNtuple*) MCFile->Get("MCValue_ntuple");
  vector<double> isotope_activity;
  float fraction;
  MCValue_ntuple->SetBranchAddress("fraction",&fraction);
  for(int j=0; j<MCValue_ntuple->GetEntries(); j++)
    {
      MCValue_ntuple->GetEntry(j);
      isotope_activity.push_back(fraction*10000.);
    }
  
  string RealDataFile = realdir + realoutfile;
  cout<<"Fitting "<<RealDataFile<<endl;
  TFile *RealFile = new TFile(RealDataFile.c_str());
  RooRealVar charge("charge","charge",0,2000);
  TTree* FullSpectrum_tree = (TTree*) RealFile->Get("lsvtree");
  RooDataSet data("data","Real Data",FullSpectrum_tree,charge);
  charge.setBins(2000,"cache");
        
  TH1F *FullSpectrum = (TH1F*) RealFile->Get("FullSpectrum");
  // RooDataHist data("data","Real Data",RooArgSet(charge),Import(*FullSpectrum));
  double integralsum = FullSpectrum->Integral();
    
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
  //  vector<RooGaussian*> mc_gauss;
  vector<RooGaussModel*> mc_gauss;
  // vector<RooProdPdf*> mc_prodpdf;
  //  vector<RooFFTConvPdf*> mc_convpdf;
  vector<RooNumConvPdf*> mc_convpdf;
  vector<RooRealVar*> mc_fraction;    
  vector<RooRealVar*> mc_rate;
  RooArgList pdflist("pdflist");
  RooArgList fraclist("fraclist");
  /*
  vector<RooRealVar*> ml_test;
  vector<RooRealVar*> sl_test;
  vector<RooLandau*> landau_test; 
  */
  vector<TCanvas*> mc_canv;
  vector<RooPlot*> mc_frame;
  vector<TCanvas*> real_canv;
  vector<RooPlot*> real_frame;

#ifdef  single_fit
  for(size_t i=iter; i<iter+1; i++)
#else
  for(size_t i=0; i<Source.size(); i++)
#endif
    {
      RooRealVar *temp_x = new RooRealVar(Form("%s_x",Source.at(i).c_str()),Form("%s x",Source.at(i).c_str()),0,NBins.at(i));
      mc_x.push_back(temp_x);      
      mc_x.back()->setBins(NBins.at(i),"cache");

      RooDataHist *temp_hist = new RooDataHist(hname.at(i),htitle.at(i),*(mc_x.back()),mcplot.at(i));    
      mc_hist.push_back(temp_hist);
      RooHistPdf *temp_pdf = new RooHistPdf(Form("%s_pdf",Source.at(i).c_str()),Form("%s pdf",Source.at(i).c_str()),*(mc_x.back()),*(mc_hist.back()),2);
      mc_pdf.push_back(temp_pdf);
      RooPlot* temp_mcframe = (*(mc_x.back())).frame();
      mc_frame.push_back(temp_mcframe);
      mc_pdf.back()->plotOn(mc_frame.back(),LineColor(linecolor.at(i)));     
      TCanvas *temp_mccanv = new TCanvas(Form("%s_mccanv",Source.at(i).c_str()),Form("%s MC Energy",Source.at(i).c_str()),600,400);
      mc_canv.push_back(temp_mccanv);
      mc_frame.back()->Draw();

      RooFormulaVar *temp_gaus_mean = new RooFormulaVar(Form("%s_gaus_mean",Source.at(i).c_str()),"@0*@1+@2",RooArgList(*(mc_x.back()),ly_mean,basel_mean));
      RooFormulaVar *temp_gaus_var = new RooFormulaVar(Form("%s_gaus_var",Source.at(i).c_str()),"@0+(1+@1)*@2*@3+@4*TMath::Power(@2*@3,2)", 
							 RooArgList(basel_var,spe_var,*(mc_x.back()),ly_mean,ly_var));    
      mc_gaus_mean.push_back(temp_gaus_mean);
      mc_gaus_var.push_back(temp_gaus_var);
      //     RooGaussian *temp_gauss = new RooGaussian(Form("%s_gauss",Source.at(i).c_str()),Form("%s Gaussian",Source.at(i).c_str()),charge,*(mc_gaus_mean.back()),*(mc_gaus_var.back()));
      RooGaussModel *temp_gauss = new RooGaussModel(Form("%s_gauss",Source.at(i).c_str()),Form("%s Gaussian",Source.at(i).c_str()),charge,*(mc_gaus_mean.back()),*(mc_gaus_var.back()));
      mc_gauss.push_back(temp_gauss);	
      
      /*
      RooProdPdf *temp_prodpdf = new RooProdPdf(Form("%s_prodpdf",Source.at(i).c_str()),Form("%s Signal & Response Prod PDF",Source.at(i).c_str()),*(mc_pdf.back()),Conditional(*(mc_gauss.back()),charge));
      mc_prodpdf.push_back(temp_prodpdf);
      pdflist.add(*(mc_prodpdf.back()));
      mc_prodpdf.back()->getObservables(*(mc_x.back()));
      */

      // RooFFTConvPdf *temp_convpdf = new RooFFTConvPdf(Form("%s_convpdf",Source.at(i).c_str()),Form("%s Signal & Response convolution",Source.at(i).c_str()),charge,*(mc_pdf.back()),*(mc_gauss.back()));
      RooNumConvPdf *temp_convpdf = new RooNumConvPdf(Form("%s_convpdf",Source.at(i).c_str()),Form("%s Signal & Response convolution",Source.at(i).c_str()),charge,*(mc_gauss.back()),*(mc_pdf.back()));
      mc_convpdf.push_back(temp_convpdf);
      // mc_convpdf.back()->setCacheObservables(*(mc_x.back())) ;
      //    mc_convpdf.back()->setBufferFraction(1.0);
      pdflist.add(*(mc_convpdf.back()));
      //RooRealVar *temp_fraction = new RooRealVar(Form("%s_frac",Source.at(i).c_str()),Form("%s fraction",Source.at(i).c_str()),isotope_activity.at(i));
      // RooRealVar temp_fraction(Form("%s_frac",Source.at(i).c_str()),Form("%s fraction",Source.at(i).c_str()),isotope_activity.at(i),0,1);
      //mc_fraction.push_back(temp_fraction);
      //fraclist.add(*(mc_fraction.back()));
      RooRealVar *temp_rate = new RooRealVar(Form("%s_rate",Source.at(i).c_str()),Form("%s Rate[Hz]",Source.at(i).c_str()),
					     integralsum*isotope_activity.at(i),0,integralsum*5);
      mc_rate.push_back(temp_rate);
      fraclist.add(*(mc_rate.back()));
   
      RooPlot* temp_frame = charge.frame();
      real_frame.push_back(temp_frame);
      //mc_prodpdf.back()->plotOn(real_frame.back(),LineColor(linecolor.at(i)));
      mc_convpdf.back()->plotOn(real_frame.back(),LineColor(linecolor.at(i)));
      TCanvas *temp_canv = new TCanvas(Form("%s_canv",Source.at(i).c_str()),Form("%s Real Energy",Source.at(i).c_str()),600,400);
      real_canv.push_back(temp_canv);
      real_frame.back()->Draw();
    }

  //    Bool_t IsRecursive = true;
  // RooAddPdf model("model","Total Fit",pdflist,fraclist,IsRecursive);	
  RooAddPdf model("model","Total Fit",pdflist,fraclist);	
  //    TH1 *model_hist = model.createHistogram();
  //  RooDataSet *datatest = model.generate(*(mc_x.back()),1000);
  // model.fitTo(*datatest);
  
  if(start_fit){
    //    model.fitTo(data,SumW2Error(kTRUE),Save(),Extended());
    //  (*(mc_prodpdf.at(0))).fitTo(data);
    RooFitResult *r = model.fitTo(data,SumW2Error(kTRUE),Save(),Extended(),NumCPU(4,0));    
    r->Print();
  }
  
  RooPlot* frame = charge.frame();
  //  RooPlot* frame = (*(mc_x.back())).frame();
  //  datatest->plotOn(frame);
  data.plotOn(frame);
  model.plotOn(frame,LineColor(linecolor.at(9)));
  //  (*(mc_prodpdf.at(0))).plotOn(frame,LineColor(linecolor.at(9)));
    
  TCanvas *c1 = new TCanvas("c1","Full Spectrum Fit",1000,600);
  gPad->SetLogy();
  frame->Draw();
  
  
  string output = dirname + "Pulse_Splitter_roofit_"+Tag+"_"+Time+".root";
  TFile outfile(output.c_str(), "RECREATE");
#ifdef  single_fit
  for(size_t i=iter; i<iter+1; i++)
#else
  for(size_t i=0; i<Source.size(); i++)
#endif
    {
      mc_canv.at(i)->Write();
      real_canv.at(i)->Write();
    }
  c1->Write();
  outfile.Write();
  outfile.Close();
    
    // }
  return 1;   
  
}
//#endif /*__CINT__*/














    /*    RooRealVar c14_rate("c14_rate","C14 Rate[Hz]",sum.at(0),0,integralsum);
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

    RooRealVar *ml = new RooRealVar(Form("ml_%s",Source.at(i).c_str()),"mean landau",5.,-20,20) ;
    RooRealVar *sl = new RooRealVar(Form("ml_%s",Source.at(i).c_str()),"sigma landau",1,0.1,10) ;
    ml_test.push_back(ml);
    sl_test.push_back(sl);
    RooLandau *lx = new RooLandau(Form("ml_%s",Source.at(i).c_str()),"lx",charge,*(ml_test.at(i)),*(sl_test.at(i))) ;
    landau_test.push_back(lx);
   
    "basel_var+(1+spe_var)*energy*ly_mean+ly_var*TMath::Power(energy*ly_mean,2)",
    RooFormulaVar temp_gaus_mean(temp_gaus_mean_name,"@0*@1+@2",RooArgList(mc_x.at(i),ly_mean,basel_mean));
    RooFormulaVar temp_gaus_var(temp_gaus_var_name,"@0+(1+@1)*@2*@3+@4*TMath::Power(@2*@3,2)", 			      
    RooArgList(basel_var,spe_var,mc_x.at(i),ly_mean,ly_var));    

    RooRealVar temp_x(temp_xname,temp_xtitle,0,NBins.at(i));
    RooDataHist *temp_hist = new RooDataHist(hname.at(i),htitle.at(i),charge,mcplot.at(i));    
    RooHistPdf *temp_pdf = new RooHistPdf(temp_pdfname,temp_pdftitle,charge,*(mc_hist.at(i)),0);

    temp_xname.Form("%s_x",Source.at(i).c_str());
    temp_xtitle.Form("%s x",Source.at(i).c_str());

    temp_pdfname.Form("%s_pdf",Source.at(i).c_str());
    temp_pdftitle.Form("%s pdf",Source.at(i).c_str());	

    temp_gaus_mean_name.Form("%s_gaus_mean",Source.at(i).c_str());
    temp_gaus_var_name.Form("%s_gaus_var",Source.at(i).c_str());
    
    RooRealVar mg("mg","mg",0) ;
    RooRealVar sg("sg","sg",2,0.1,10) ;
    RooGaussian *temp_gauss = new RooGaussian(temp_gauss_name,temp_gauss_title,charge,mg,sg);

    temp_gauss_name.Form("%s_gauss",Source.at(i).c_str());
    temp_gauss_title.Form("%s Gaussian",Source.at(i).c_str());

    temp_convpdfname.Form("%s_convpdf",Source.at(i).c_str());
    temp_convpdftitle.Form("%s Signal & Response convolution",Source.at(i).c_str());	
    RooFFTConvPdf *temp_convpdf = new RooFFTConvPdf(temp_convpdfname,temp_convpdftitle,charge,*(landau_test.at(i)),*(mc_gauss.at(i)));
    temp_prodpdfname.Form("%s_prodpdf",Source.at(i).c_str());
    temp_prodpdftitle.Form("%s Signal & Response Prod PDF",Source.at(i).c_str());	
    RooProdPdf *temp_prodpdf = new RooProdPdf(temp_prodpdfname,temp_prodpdftitle,charge,*(landau_test.at(i)),*(mc_gauss.at(i)));
    TString temp_xname,temp_xtitle,temp_pdfname,temp_pdftitle,temp_gaus_mean_name,temp_gaus_var_name,
    temp_gauss_name,temp_gauss_title,temp_convpdfname,temp_convpdftitle,temp_prodpdfname,temp_prodpdftitle;

    RooPlot* frame = charge.frame();
    data.plotOn(frame);
    TCanvas *c1 = new TCanvas("c1","Full Spectrum Fit",1000,600);
    frame->Draw();

    //#ifndef __CINT__
    int main(int argc, char** argv){
    theApp = new TRint("App",&argc,argv,NULL,0);
    bool start_fit;
    int tag;
    if(theApp->Argc() == 2)
    {
    tag = theApp->Argc();
    start_fit = atoi(theApp->Argv(1))>0 ? true : false;
    }
    else{
    cout<<"Usage: ./reconMain start_fit "<<endl;
    return 0;
    }



    */
