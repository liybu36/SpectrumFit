#include "vetopulsesplit.hh"
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "TH1F.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TRint.h"
#include "TMath.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TROOT.h"
#include "TCut.h"
#include "TLegend.h"
#include "TPad.h"
#include "TNtuple.h"

using namespace std;
ClassImp(vetopulsesplit);

double vetopulsesplit::Isotope_Sum(int isotopes, vector<double> isotope_activity)
{
  double activity_sum=0;
  for(int i=0; i<isotopes; i++)    
    activity_sum += isotope_activity.at(i);
  return activity_sum;
}

void vetopulsesplit::Isotope_Fraction(int isotopes, vector<double> isotope_activity, double *isotope_fraction)
{
  double sum = Isotope_Sum(isotopes,isotope_activity);
  for(int i=0; i<isotopes; i++)   
    {
      isotope_fraction[i] = isotope_activity.at(i)*1.0/sum; 
    }
}

void vetopulsesplit::Calculate_Convolution(TH1F *energy_spectrum1, TH1F *energy_spectrum2, TH1F *convolution)
{
  int np1 = energy_spectrum1->GetNbinsX();
  int np2 = energy_spectrum2->GetNbinsX();
  int nconv = np1 + np2;
  double binsize;
  for (int n = 1; n <= nconv; n++) {
    double sum = 0;
    for (int m = 1; m <= n; m++) {
      if (m > np1  || (n - m) > np2 ) continue;
      sum += (energy_spectrum1->GetBinContent(m)) * (energy_spectrum2->GetBinContent(n-m+1));
    }
    convolution->SetBinContent(n,sum);
   }
  double integral = convolution->Integral();
  cout<<convolution->GetName()<<"\t integral= "<<integral<<endl;
  convolution->Scale(1./integral);
}

void vetopulsesplit::Calculate_Norm (TH1F *energy_spectrum,  vector<double> &convolution_spectrum)
{
  int np = energy_spectrum->GetNbinsX();
  double integral = energy_spectrum->Integral();
  cout<<"integral= "<<integral<<endl;
  energy_spectrum->Scale(1./integral);
  for(int i=1; i<=np; i++)
    convolution_spectrum.push_back(energy_spectrum->GetBinContent(i));
}

double vetopulsesplit::scint_e_quenching (double e_keV, double* params) {

  double kB =0;// params[1];


  double birks_value[7][6] = {{0.49041,0.21600,0.10725,-1.46804e-3,0.22468,0.09544}, //0.006
			      {0.42149,0.18731,0.09624,-0.52346e-3,0.19977,0.08311}, //0.008
			      {0.36851,0.16398,0.08816,-0.45562e-3,0.17582,0.07431}, //0.010
			      {0.32903,0.14404,0.08059,-0.04536e-3,0.15623,0.06611}, //0.012
			      {0.29668,0.12872,0.07477,0.15992e-3,0.14091,0.05978}, //0.014
			      {0.27020,0.11663,0.06985,0.34938e-3,0.12933,0.05453}, //0.016
			      {0.24808,0.10646,0.06576,0.46844e-3,0.11938,0.04998}  //0.018
  };
  

  double A1 = -0.6292-0.2181*log(kB);
  double A2 = -0.3057-0.1024*log(kB);
  double A3 = -0.0673-0.03353*log(kB);
  double A4 = 0.009876+0.002276*log(kB);
  double A5 = -0.2814-0.09964*log(kB);
  double A6 = -0.09962-0.0376*log(kB);
  
  return ( (A1 + A2*log(e_keV) + A3*log(e_keV)*log(e_keV) +
	    A4*log(e_keV)*log(e_keV)*log(e_keV))/(1 + A5*log(e_keV) 
	   + A6*log(e_keV)*log(e_keV) + A4*log(e_keV)*log(e_keV)*log(e_keV)) );
}

double vetopulsesplit::response_function (double q, double energy, double* params, int choice) {
  double ly_mean = params[1];
  double spe_var = params[2];
  double ly_var = params[3];
  double baseline_mean = params[5];
  double baseline_var = params[6];
  double threshold = params[7];
  double pe;
  double gaus_mean;

  pe = energy*ly_mean;
  //  if (energy == 0) pe = energy*ly_mean;
  //pe = energy*ly_mean*scint_e_quenching(energy, params);

  if (choice == 1) gaus_mean = pe + (baseline_mean + threshold);
  else gaus_mean = pe + baseline_mean;

  double gaus_var = baseline_var + (1 + spe_var)*pe + ly_var*pe*pe;
  //  double gaus_var = baseline_var*(1+spe_var)*gaus_mean+ly_var*TMath::Power(gaus_mean,2);
  double gaus_var_inv = 1.0 / gaus_var;
  double arg = -0.5*(q - gaus_mean)*(q - gaus_mean)*gaus_var_inv;

  return 0.3989422804 * sqrt(gaus_var_inv) * exp(arg);
}

double vetopulsesplit::C14Fit (double* x, double* params){
  int nBins = EnergyMC[0]->GetNbinsX();
  double delta = EnergyMC[0]->GetBinWidth(1);
  int choice = 1;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[0]->GetBinCenter(i);
      spectrum = EnergyMC[0]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double C14Rate = params[9];
  //  return result*(C14Rate)*delta;
  //  return result*(C14Rate*(1./6103.64)*1.e-7*177353.)*params[0]*delta;

  return result*(C14Rate*(2100.e-9)/449177)*params[0]*delta; //fieldon
  //  return result*(C14Rate*(2100.e-9)/253431)*params[0]*delta; //fieldoff
}

double vetopulsesplit::C14Fit_1st (double* x, double* params){
  int nBins = C14_EnergyMC_1st->GetNbinsX();
  double delta = C14_EnergyMC_1st->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  double window = params[8];
  for (int i=1; i<=nBins; i++)
    {
      energy = C14_EnergyMC_1st->GetBinCenter(i);
      spectrum = C14_EnergyMC_1st->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double stRate = pow(params[9],2)*window*exp(-params[9]*window);
  //  return result*(stRate)*delta;
  return result*(stRate)*params[0]*delta;
}

double vetopulsesplit::C14Fit_2nd (double* x, double* params){
  int nBins = C14_EnergyMC_2nd->GetNbinsX();
  double delta = C14_EnergyMC_2nd->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  double window = params[8];
  for (int i=1; i<=nBins; i++)
    {
      energy = C14_EnergyMC_2nd->GetBinCenter(i);
      spectrum = C14_EnergyMC_2nd->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double ndRate = pow(params[9],2)*window*exp(-params[9]*window);
  //  return result*(ndRate)*delta;
  return result*(ndRate)*params[0]*delta;
}

double vetopulsesplit::Co60Fit (double* x, double* params){
  int nBins = EnergyMC[1]->GetNbinsX();
  double delta = EnergyMC[1]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[1]->GetBinCenter(i);
      spectrum = EnergyMC[1]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double Co60Rate = params[10];
  //  return result*(Co60Rate)*delta;
  return result*(Co60Rate)*params[0]*delta;
}

double vetopulsesplit::Co57Fit (double* x, double* params){
  int nBins = EnergyMC[2]->GetNbinsX();
  double delta = EnergyMC[2]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[2]->GetBinCenter(i);
      spectrum = EnergyMC[2]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double Co57Rate = params[11];
  //  return result*(Co57Rate)*delta;
  return result*(Co57Rate)*params[0]*delta;
}

double vetopulsesplit::K40Fit (double* x, double* params){
  int nBins = EnergyMC[3]->GetNbinsX();
  double delta = EnergyMC[3]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[3]->GetBinCenter(i);
      spectrum = EnergyMC[3]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double K40Rate = params[12];
  //  return result*(K40Rate)*delta;
  return result*(K40Rate)*params[0]*delta;
}

double vetopulsesplit::Tl208Fit (double* x, double* params){
  int nBins = EnergyMC[4]->GetNbinsX();
  double delta = EnergyMC[4]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[4]->GetBinCenter(i);
      spectrum = EnergyMC[4]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double Tl208Rate = params[13];
  //  return result*(Tl208Rate)*delta;
  return result*(Tl208Rate)*params[0]*delta;
}

double vetopulsesplit::Th232Fit (double* x, double* params){
  int nBins = EnergyMC[4]->GetNbinsX();
  double delta = EnergyMC[4]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[4]->GetBinCenter(i);
      spectrum = EnergyMC[4]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double Th232Rate = params[13];
  //  return result*(Th232Rate)*delta;
  return result*(Th232Rate)*params[0]*delta;
}

double vetopulsesplit::Th232LowerFit (double* x, double* params){
  int nBins = EnergyMC[5]->GetNbinsX();
  double delta = EnergyMC[5]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[5]->GetBinCenter(i);
      spectrum = EnergyMC[5]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double Th232LowerRate = params[14];
  return result*(Th232LowerRate)*params[0]*delta;
}

double vetopulsesplit::U235Fit (double* x, double* params){
  int nBins = EnergyMC[6]->GetNbinsX();
  double delta = EnergyMC[6]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[6]->GetBinCenter(i);
      spectrum = EnergyMC[6]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }

  double U235Rate = params[17]/21.5;
  return result*(U235Rate)*params[0]*delta;
}

double vetopulsesplit::U235LowerFit (double* x, double* params){
  int nBins = EnergyMC[7]->GetNbinsX();
  double delta = EnergyMC[7]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[7]->GetBinCenter(i);
      spectrum = EnergyMC[7]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double U235LowerRate = params[16];
  return result*(U235LowerRate)*params[0]*delta;
}

double vetopulsesplit::U238UpperFit (double* x, double* params){
  int nBins = EnergyMC[8]->GetNbinsX();
  double delta = EnergyMC[8]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[8]->GetBinCenter(i);
      spectrum = EnergyMC[8]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double U238UpperRate = params[17];
  //  return result*(U238UpperRate)*delta;
  return result*(U238UpperRate)*params[0]*delta;
}

double vetopulsesplit::U238LowerFit (double* x, double* params){
  int nBins = EnergyMC[9]->GetNbinsX();
  double delta = EnergyMC[9]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[9]->GetBinCenter(i);
      spectrum = EnergyMC[9]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double U238LowerRate = params[18];
  //  return result*(U238LowerRate)*delta;
  return result*(U238LowerRate)*params[0]*delta;
}

double vetopulsesplit::U238Fit (double* x, double* params){
  int nBins = EnergyMC[6]->GetNbinsX();
  double delta = EnergyMC[6]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[6]->GetBinCenter(i);
      spectrum = EnergyMC[6]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double U238Rate = params[15];
  return result*(U238Rate)*params[0]*delta;
}

double vetopulsesplit::Rn222Fit (double* x, double* params){
  int nBins = EnergyMC[6]->GetNbinsX();
  double delta = EnergyMC[6]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[6]->GetBinCenter(i);
      spectrum = EnergyMC[6]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double Rn222Rate = params[15];
  return result*(Rn222Rate)*params[0]*delta;
}

double vetopulsesplit::K42Fit (double* x, double* params){
  int nBins = EnergyMC[8]->GetNbinsX();
  double delta = EnergyMC[8]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[8]->GetBinCenter(i);
      spectrum = EnergyMC[8]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double K42Rate = params[17];
  return result*(K42Rate)*params[0]*delta;
}

double vetopulsesplit::TotalFit(Double_t *x, Double_t *params)
{
  double C14_value = C14Fit(x,params);
  //double C14_1st_value = C14Fit_1st(x,params);
  //  double C14_2nd_value = C14Fit_2nd(x,params);
  double Co60_value = Co60Fit(x,params);
  double Co57_value = Co57Fit(x,params);
  double K40_value = K40Fit(x,params);
  //  double Tl208_value = Tl208Fit(x,params);
  double Th232_value = Th232Fit(x,params);
  double Th232Lower_value = Th232LowerFit(x,params);
  double U235_value = U235Fit(x,params); 
  double U235Lower_value = U235LowerFit(x,params); 
  double U238Upper_value = U238UpperFit(x,params);
  double U238Lower_value = U238LowerFit(x,params);
  // double U238_value = U238Fit(x,params);
  //  double Rn222_value = Rn222Fit(x,params);
  // double K42_value = K42Fit(x,params);
    
  // return C14_value+Co60_value+Co57_value+K40_value+Th232_value+U238Upper_value+U238Lower_value+U235_value+K42_value+params[4];
  return C14_value+Co60_value+Co57_value+K40_value+Th232_value+U238Upper_value+U238Lower_value+U235_value+Th232Lower_value+U235Lower_value+params[4];
  //  return C14_value + Co60_value+Co57_value+K40_value+Th232_value+U238_value+K42_value+U235_value+params[4];
  //return C14_value + Co60_value+Co57_value+K40_value+Th232_value+Rn222_value+K42_value+U235_value+params[4];  
}

void vetopulsesplit::SetIsotopes()
{
  Source.push_back("C14");
  Source.push_back("Co60");
  Source.push_back("Co57");
  Source.push_back("K40");
  Source.push_back("Th232");
  Source.push_back("Th232Lower");
  Source.push_back("U235");  
  Source.push_back("U235Lower");  
  Source.push_back("U238");
  Source.push_back("U238Lower");
  //Source.push_back("Rn222");
  //Source.push_back("K42");
  //Source.push_back("Ni65");  
  //Source.push_back("Fe59");  
  
  PDG.push_back(1000060140);
  PDG.push_back(1000270600);
  PDG.push_back(1000270570);
  PDG.push_back(1000190400);
  PDG.push_back(1000902320);
  PDG.push_back(1000862200);
  PDG.push_back(1000922350);
  PDG.push_back(1000862190);
  PDG.push_back(1000922380);
  PDG.push_back(1000862220);
  //PDG.push_back(1000190420);
  //PDG.push_back(1000280650);
  //PDG.push_back(1000260590);
  
  Normalization.resize(Source.size());
  fraction.resize(Source.size());

  NBins.push_back(160);
  NBins.push_back(2600);
  NBins.push_back(800);
  NBins.push_back(1600);
  NBins.push_back(2000);
  NBins.push_back(3500);
  NBins.push_back(900);
  NBins.push_back(1200);
  NBins.push_back(1800);
  NBins.push_back(3000);
  //NBins.push_back(3500);
  //NBins.push_back(3000);
  //NBins.push_back(1600);
  //NBins.push_back(1400);

  //TPC energy range
  TBins.push_back(160);
  TBins.push_back(1600);
  TBins.push_back(500);
  TBins.push_back(1000);
  TBins.push_back(1300);
  TBins.push_back(2500);
  TBins.push_back(500);
  TBins.push_back(800);
  TBins.push_back(1200);
  TBins.push_back(2000);

}

bool vetopulsesplit::multicut(float height,float multiplicity, float charge){
  return height/multiplicity < (2.563e7 + TMath::Sqrt(1.574e14+1.390e12*(charge-14.40)*(charge-14.40)));
}

void vetopulsesplit::Readdatafile(TChain *t, int startfile, int endfile)
{
  //  string dirname="/darkside/users/hqian/pulse_splitter/";
  string dirname=GetRealinputdir();
  string filename;
  string middle="ODRun000";
  string last=".root";
  stringstream oss;
  ifstream NameCheck;

  for(int i=startfile; i<= endfile; i++)
    {
      oss<<i;
      filename=dirname+middle+oss.str()+last;
      NameCheck.open(filename.c_str());
      if(!NameCheck.good())
	{
	  oss.str("");
	  NameCheck.close();
	  continue;
	}
      else{
	t->Add(filename.c_str());
	cout<<"Processing Data file: "<<filename<<endl;
	oss.str("");
      }
      NameCheck.close();
    }
}

void vetopulsesplit::DST_Readdatafile(TChain *t, int startfile, int endfile)
{
  bool fieldon =true;
  if(fieldon)
    {
      //run 11856 to 11928
      string dirname="/darkside/data/UAr_DSTs/";
      //  string dirname=GetRealinputdir();
      string middle="DST_Run";
      string last=".root";
      for(int i=startfile; i<=endfile; i++)
	{
	  TString filename = Form("%s%06d%s",middle.c_str(),i,last.c_str());
	  filename.Prepend(dirname.c_str());
	  ifstream NameCheck;
	  NameCheck.open(filename.Data());
	  if(!NameCheck.good())
	    continue;
	  else{
	    TFile *f = new TFile(filename);
	    if(f->IsZombie())
	      continue;
	    else{
	      t->Add(filename);
	      cout<<"Processing Data file: "<<filename<<endl;
	    }
	  }
	}
    }
  else{
    string dirname="/darkside/data/UAr_DSTs/nullfield/";
    string filename = dirname + "DST_Run011764_11822.root";
    TFile *f = new TFile(filename.c_str());
    if(!f->IsZombie())
      {
	t->Add(filename.c_str());
	cout<<"Processing Data file: "<<filename<<endl;
      }
  }
}

void vetopulsesplit::Recon_Readdatafile(TChain *t, int startfile, int endfile,int k,string last)
{
  string dirname=GetMCinputdir();
  string middle="out"+Source.at(k);
  //  string last="_clustered.root";
  for(int i=startfile; i<=endfile; i++)
    {
      TString filename;
      if(i==0) filename.Form("%s%s",middle.c_str(),last.c_str());
      else
	filename.Form("%s_v%d%s",middle.c_str(),i,last.c_str());
      filename.Prepend(dirname.c_str());
      ifstream NameCheck;
      NameCheck.open(filename.Data());
      if(!NameCheck.good())
	continue;
      else{
	TFile *f = new TFile(filename);
	if(f->IsZombie())
	  continue;
	else{
	  t->Add(filename);
	  cout<<"Processing Data file: "<<filename<<endl;
	}
      }
    }
}

//double vetopulsesplit::Calculate_MCData(TTree *dstree, TH1F *Energy_Spectrum)
void vetopulsesplit::Calculate_MCData(int k, TH1F *Energy_Spectrum)
{
  string MCinputdir=GetMCinputdir();  
  string MCinputfile = MCinputdir+"out" + Source.at(k)+ "_v4_clustered.root";
  TFile *file = new TFile(MCinputfile.c_str());
  TTree *dstree = (TTree *) file->Get("Recon");
  if(!file->IsOpen())
    cout<<"Error: Cannot Open MC Data File =>"<<MCinputfile<<endl;       
  else  cout<<"Processing MC Data file: "<<MCinputfile<<endl;  

  Int_t NEntries = dstree->GetEntries();
  cout<<"NEntries= "<<NEntries<<endl;
  Double_t  sum;

  vector<TString> *volume;
  vector<double>  *et;
  vector<double>  *ex;
  vector<double>  *ey;
  vector<double>  *ez;
  vector<double>  *edep;
  vector<double>  *eqch;
  vector<double>  *quenchingfactor;
  volume = 0;
  et = 0;
  ex = 0;
  ey = 0;
  ez = 0;
  edep = 0;
  eqch = 0;
  quenchingfactor = 0;
  dstree->SetBranchAddress("volume",&volume);
  dstree->SetBranchAddress("et",&et);
  dstree->SetBranchAddress("ex",&ex);
  dstree->SetBranchAddress("ey",&ey);
  dstree->SetBranchAddress("ez",&ez);
  dstree->SetBranchAddress("edep",&edep);
  dstree->SetBranchAddress("eqch",&eqch);
  dstree->SetBranchAddress("quenchingfactor",&quenchingfactor);

  int counts=0;
  for(Int_t i=0; i<NEntries; i++)
    { dstree->GetEntry(i);
      bool capture = false;
      for(Int_t j=0; j<et->size(); j++)
	{
	  if(volume->at(j)=="p_scint")
	    { capture = true;
	      Energy_Spectrum->Fill(eqch->at(j));
	      EdepMC.at(k)->Fill(edep->at(j));
	    }
	}
      if(capture) ++counts;  
    }
  Int_t NBins = Energy_Spectrum->GetNbinsX();
  Double_t integralsum = Energy_Spectrum->Integral();
  cout<<"NBins= "<<NBins<<"\t integralsum="<<integralsum<<"\t counts="<<counts<<endl; 
  //  Energy_Spectrum->Scale(1./integralsum);
  //  Energy_Spectrum->Scale(1./NEntries);
  Energy_Spectrum->Scale(1./Normalization.at(k));
  EdepMC.at(k)->Scale(1./Normalization.at(k));
  //  return (1.0*counts/NEntries);
  fraction.at(k) = 1.0*counts/NEntries;
}

bool Compare_Time(double a, double b)
{
  return a<b;
}

//void vetopulsesplit::Calculate_MCCoincidenceData(TTree *dstree, TH1F *Energy_Spectrum, TNtuple *veto_ntuple)
//double vetopulsesplit::Calculate_MCCoincidenceData(TTree *dstree, TH1F *Energy_Spectrum)
void vetopulsesplit::Calculate_MCCoincidenceData(int k, TH1F *Energy_Spectrum)
{
  string last="_clustered.root";
  TChain *dstree = new TChain("Recon");
  Recon_Readdatafile(dstree,4,4,k,last);

  Int_t NEntries = dstree->GetEntries();
  cout<<"NEntries= "<<NEntries<<endl;
  Double_t  sum;

  //  Int_t           pdg;
  //  Bool_t          event_broken;
  vector<TString> *volume;
  vector<double>  *et;
  vector<double>  *ex;
  vector<double>  *ey;
  vector<double>  *ez;
  vector<double>  *edep;
  vector<double>  *eqch;
  vector<double>  *quenchingfactor;
  volume = 0;
  et = 0;
  ex = 0;
  ey = 0;
  ez = 0;
  edep = 0;
  eqch = 0;
  quenchingfactor = 0;
  //  dstree->SetBranchAddress("event_pdg", &pdg);
  //  dstree->SetBranchAddress("event_broken", &event_broken);
  dstree->SetBranchAddress("volume",&volume);
  dstree->SetBranchAddress("et",&et);
  dstree->SetBranchAddress("ex",&ex);
  dstree->SetBranchAddress("ey",&ey);
  dstree->SetBranchAddress("ez",&ez);
  dstree->SetBranchAddress("edep",&edep);
  dstree->SetBranchAddress("eqch",&eqch);
  dstree->SetBranchAddress("quenchingfactor",&quenchingfactor);
  
  double tpc_low_threshold=0.5;
  //  double tpc_high_threshold=600;
  double prompt_time = -50.; //ns
  double delay_time = 50.; //ns
  //  std::map<Int_t,Int_t> pdg_counts; 
  int counts=0;
  for(Int_t i=0; i<NEntries; i++)
    { dstree->GetEntry(i);
      /*      if(pdg_counts.count(pdg)==0)	      
	pdg_counts.insert( std::pair<Int_t,Int_t>(pdg,1) );
      else{
	++(pdg_counts.find(pdg)->second);	
      */
      double tpc_total=0;
      vector<double> tpc_trigger_time;
      bool capture = false;     
      for(Int_t j=0; j<et->size(); j++)
	{
	  if(volume->at(j)=="p_active" && eqch->at(j)>tpc_low_threshold)// && eqch->at(j)<tpc_high_threshold)
	    {
	      tpc_trigger_time.push_back(et->at(j)*1.e+9);	   
	      TPC_EnergyMC.at(k)->Fill(eqch->at(j));
	      TPC_EdepMC.at(k)->Fill(edep->at(j));
	      tpc_total += eqch->at(j);
	    }
	}
      //#define Coincidence 
#ifdef Coincidence
      if(tpc_trigger_time.size()>1)
	std::sort(tpc_trigger_time.begin(),tpc_trigger_time.end(),Compare_Time);
      if(tpc_trigger_time.size())
	{ 
	  for(Int_t j=0; j<et->size(); j++)
	    {
	      if(volume->at(j)=="p_scint")
		{ double gps = et->at(j)*1.e+9 - tpc_trigger_time.front();
		  // if(gps>prompt_time && gps<delay_time)
		  if(gps>-2000. && gps<100.)
		    { capture = true;
		      Energy_Spectrum->Fill(eqch->at(j));		 
		      EdepMC.at(k)->Fill(edep->at(j));
		      TPC_Veto_EnergyMC.at(k)->Fill(tpc_total,eqch->at(j));
		    }  	}   }	}
#else
      for(Int_t j=0; j<et->size(); j++)
	{
	  if(volume->at(j)=="p_scint")
	    {
	      capture = true;
	      Energy_Spectrum->Fill(eqch->at(j));		 
	      EdepMC.at(k)->Fill(edep->at(j));
	    }
	}
#endif
      if(capture) ++counts;
    }	  
  Energy_Spectrum->Rebin(2);
  Int_t NBins = Energy_Spectrum->GetNbinsX();
  Double_t integralsum = Energy_Spectrum->Integral();
  cout<<"NBins= "<<NBins<<"\t integralsum="<<integralsum<<"\t counts="<<counts<<endl; 
  //  Energy_Spectrum->Scale(1./integralsum);
  //  Energy_Spectrum->Scale(1./NEntries);
  /*
  for (std::map<Int_t, Int_t>::iterator it=pdg_counts.begin(); it!=pdg_counts.end(); ++it)
    cout << it->first << " => " << it->second << '\n';
  Normalization.at(k) = pdg_counts.find(PDG.at(k))->second;
  cout<<"Normalization= "<<Normalization.at(k)<<endl;    
  */
  Energy_Spectrum->Scale(1./Normalization.at(k));
  EdepMC.at(k)->Scale(1./Normalization.at(k));
  //  return (1.0*counts/NEntries);
  fraction.at(k) = 1.0*counts/NEntries;
}

void vetopulsesplit::Calculate_U238MCCoincidenceData(int k, TH1F *Energy_Spectrum, TH1F *Next)
{
  string last="_clustered.root";
  TChain *dstree = new TChain("Recon");
  Recon_Readdatafile(dstree,4,4,k,last);

  Int_t NEntries = dstree->GetEntries();
  cout<<"NEntries= "<<NEntries<<endl;
  Double_t  sum;

  Int_t           pdg;
  Bool_t          event_broken;
  vector<TString> *volume;
  vector<double>  *et;
  vector<double>  *ex;
  vector<double>  *ey;
  vector<double>  *ez;
  vector<double>  *edep;
  vector<double>  *eqch;
  vector<double>  *quenchingfactor;
  volume = 0;
  et = 0;
  ex = 0;
  ey = 0;
  ez = 0;
  edep = 0;
  eqch = 0;
  quenchingfactor = 0;
  dstree->SetBranchAddress("event_pdg", &pdg);
  dstree->SetBranchAddress("event_broken", &event_broken);
  dstree->SetBranchAddress("volume",&volume);
  dstree->SetBranchAddress("et",&et);
  dstree->SetBranchAddress("ex",&ex);
  dstree->SetBranchAddress("ey",&ey);
  dstree->SetBranchAddress("ez",&ez);
  dstree->SetBranchAddress("edep",&edep);
  dstree->SetBranchAddress("eqch",&eqch);
  dstree->SetBranchAddress("quenchingfactor",&quenchingfactor);
  
  double tpc_low_threshold=0.5;
  //  double tpc_high_threshold=600;
  double prompt_time = -50.; //ns
  double delay_time = 50.; //ns

  for(Int_t i=0; i<NEntries; i++)
    { dstree->GetEntry(i);
      vector<double> tpc_trigger_time;
      for(Int_t j=0; j<et->size(); j++)
	{
	  if(volume->at(j)=="p_active" && eqch->at(j)>tpc_low_threshold)// && eqch->at(j)<tpc_high_threshold)
	    {
	      tpc_trigger_time.push_back(et->at(j)*1.e+9);	   
	      if(!event_broken)
		{
		  TPC_EnergyMC.at(k)->Fill(eqch->at(j));
		  TPC_EdepMC.at(k)->Fill(edep->at(j));
		}
	      else
		{
		  TPC_EnergyMC.at(k+1)->Fill(eqch->at(j));
		  TPC_EdepMC.at(k+1)->Fill(edep->at(j));
		}
	    }
	}
      //#define U238_Coincidence
#ifdef U238_Coincidence
      if(tpc_trigger_time.size()>1)
	std::sort(tpc_trigger_time.begin(),tpc_trigger_time.end(),Compare_Time);
      if(tpc_trigger_time.size())
	{ 
	  for(Int_t j=0; j<et->size(); j++)
	    {
	      if(volume->at(j)=="p_scint")
		{ double gps = et->at(j)*1.e+9 - tpc_trigger_time.front();
		  if(gps>-2000. && gps<100.)
		    //  if(gps>prompt_time && gps<delay_time)
		    { 
		      if(!event_broken)
			{
			  Energy_Spectrum->Fill(eqch->at(j));
			  EdepMC.at(k)->Fill(edep->at(j));
			}
		      else{
			Next->Fill(eqch->at(j));
			EdepMC.at(k+1)->Fill(edep->at(j));
		      }	 } 
		} } }
#else
      for(Int_t j=0; j<et->size(); j++)
	{
	  if(volume->at(j)=="p_scint")
	    {	
	      if(!event_broken){
		EdepMC.at(k)->Fill(edep->at(j));
		Energy_Spectrum->Fill(eqch->at(j));
	      }
	      else{
		Next->Fill(eqch->at(j));	      
		EdepMC.at(k+1)->Fill(edep->at(j));
	      }
	    }
	}
#endif
    }

  Energy_Spectrum->Rebin(2);
  Next->Rebin(2);
  Int_t NBins = Energy_Spectrum->GetNbinsX();
  Double_t integralsum[2];
  integralsum[0] = Energy_Spectrum->Integral();
  integralsum[1] = Next->Integral();
  for(int j=0; j<2; j++){
    fraction.at(k+j) = 1.0*integralsum[j]/NEntries;
  }
 
  Energy_Spectrum->Scale(1./Normalization.at(k));
  Next->Scale(1./Normalization.at(k+1));
  EdepMC.at(k)->Scale(1./Normalization.at(k));
  EdepMC.at(k+1)->Scale(1./Normalization.at(k+1));  
}

void vetopulsesplit::Dstree_Analysis(int k)
{ 
  string last=".root";
  TChain *dstree = new TChain("dstree");
  Recon_Readdatafile(dstree,0,6,k,last);

  Int_t NEntries = dstree->GetEntries();
  cout<<"NEntries= "<<NEntries<<endl;
  
  Int_t           pdg;
  dstree->SetBranchAddress("pdg", &pdg);
  std::map<Int_t,Int_t> pdg_counts; 
  for(int i=0; i<NEntries; i++)
    {
      dstree->GetEntry(i);
      if(pdg_counts.count(pdg)==0)	      
	pdg_counts.insert( std::pair<Int_t,Int_t>(pdg,1) );
      else{
	++(pdg_counts.find(pdg)->second);	
      }
    }
  for (std::map<Int_t, Int_t>::iterator it=pdg_counts.begin(); it!=pdg_counts.end(); ++it)
    cout << it->first << " => " << it->second << '\n';
  Normalization.at(k) = pdg_counts.find(PDG.at(k))->second;
  cout<<"Normalization= "<<Normalization.at(k)<<endl;    
  
  //  if(k==6)
  if(k>=4)
    {
      Normalization.at(k+1) = pdg_counts.find(PDG.at(k+1))->second;
      cout<<"Normalization= "<<Normalization.at(k+1)<<endl;       
    } 
}

void vetopulsesplit::Recon_DataAnalysis()
{
  //  SetIsotopes();
  string MCinputdir=GetMCinputdir();
  string temp = "_EnergyMC";
  for(size_t i=0; i<Source.size(); i++)  
  {
    TString hname=Form("%s%s",Source[i].c_str(),temp.c_str());
    TString htitle=Form("%s MC Energy",Source[i].c_str());
    EnergyMC.push_back(new TH1F(hname,htitle,NBins[i],0,NBins[i]*1.0));
    EnergyMC.at(i)->GetXaxis()->SetTitle("Quenching Energy [keVee]");
    Histlist.Add(EnergyMC[i]);      
    EdepMC.push_back(new TH1F(Form("%s_EdepMC",Source[i].c_str()),"; Energy [keV]",NBins[i],0,NBins[i]*1.0));
    Histlist.Add(EdepMC[i]);      
    
    TString tpcname=Form("%s_TPC%s",Source[i].c_str(),temp.c_str());
    TString tpctitle=Form("%s TPC MC Energy",Source[i].c_str());
    TPC_EnergyMC.push_back(new TH1F(tpcname,tpctitle,TBins[i],0,TBins[i]*1.0));
    Histlist.Add(TPC_EnergyMC[i]);      
    TPC_EdepMC.push_back(new TH1F(Form("%s_TPC_EdepMC",Source[i].c_str()),"; Energy [keV]",TBins[i],0,TBins[i]*1.0));
    Histlist.Add(TPC_EdepMC[i]);      

    TPC_Veto_EnergyMC.push_back(new TH2F(Form("%s_TPC_Veto%s",Source[i].c_str(),temp.c_str()),Form("%s TPC Veto MC Energy",Source[i].c_str()),NBins[i],0,NBins[i]*1.0,NBins[i],0,NBins[i]*1.0));
    TPC_Veto_EnergyMC.at(i)->GetXaxis()->SetTitle("TPC Energy [keV]");
    TPC_Veto_EnergyMC.at(i)->GetYaxis()->SetTitle("Veto Energy [keV]");
    Histlist.Add(TPC_Veto_EnergyMC[i]);      
  }
  MCValue_ntuple = new TNtuple("MCValue_ntuple","MC Stats Value","fraction:Normalization");
  Histlist.Add(MCValue_ntuple);
  for(size_t i=0; i<Source.size(); i++)    
    {
      if(i==0)
      //if(i==0 || i==1 || i==3)
	{
	  Dstree_Analysis(i);	
	  Calculate_MCData(i,EnergyMC[i]);
	}
      else
	{	  
	  // if(i==6)
	  if(i>=4)
	  {
	      Dstree_Analysis(i);
	      Calculate_U238MCCoincidenceData(i,EnergyMC[i],EnergyMC[i+1]);    
	      ++i;
	    }
	  else{	
	    Dstree_Analysis(i);
	    Calculate_MCCoincidenceData(i,EnergyMC[i]);    	 	  
	  }
	} 
    } 
  for(size_t i=0; i<Source.size(); i++)    
    MCValue_ntuple->Fill(fraction.at(i),Normalization.at(i));  
 
  string outdir=GetMCinputdir();
  //  string output = outdir +"PulseSplitMCEnergy"+label+".root";
  string output = outdir + GetMCFile();
  TFile outfile(output.c_str(), "RECREATE");
  Histlist.Write();                       
  outfile.Write();
  outfile.Close();
  cout<<"Successfully save MC Data."<<endl;
}  

void vetopulsesplit::ODTree_DataAnalysis(int start, int end)
{
  TChain *odtree =new TChain("odtree");
  Readdatafile(odtree,start,end);
  Int_t nEntries = odtree->GetEntries();
  cout<<"nEntries= "<<nEntries<<endl;

  string outdir=GetRealinputdir();
#ifdef MoreBins
  string output = outdir +"PulseSplitRealEnergy_MoreBins.root";
  int Bins = 2000;
#else
  string output = outdir +"PulseSplitRealEnergy.root";
  int Bins = 1000;
#endif

  TH1F* FullSpectrum = new TH1F("FullSpectrum","full spectrum;charge [PE];Counts",Bins,0,2000);
  TCut TimingCut = "lsv_cluster_start_ns.fArray>3770 && lsv_cluster_start_ns.fArray<3786";
  TCut MultiCut = "lsv_cluster_height.fArray/lsv_cluster_max_multiplicity.fArray < (2.536e+7 + TMath::Sqrt(1.574e14+1.390e12*(lsv_cluster_fixed_width_charge.fArray-14.40)**2))";
  TCanvas *c1 = new TCanvas("c1", "Full Spectrum with Fit", 1200, 600);
  c1->SetLogy();
  c1->cd();
  odtree->Draw("lsv_cluster_fixed_width_charge.fArray>>FullSpectrum", TimingCut && MultiCut && "lsv_cluster_fixed_width_charge.fArray>0");
#define MOreBins
#ifndef MoreBins
  FullSpectrum->Sumw2();
  FullSpectrum->Scale(1.0/nEntries);
#endif

  TFile outfile(output.c_str(), "RECREATE");
  FullSpectrum->Write();
  c1->Write();
  outfile.Write();
  outfile.Close();
  cout<<"Successfully save Real Data."<<endl;

}
































/*
void vetopulsesplit::ODTree_Branch(TChain *fChain)
{
  lsv_cluster_fixed_width_charge = 0;
  lsv_cluster_sat_corr_charge = 0;
  lsv_cluster_start_ns = 0;
  lsv_cluster_width_ns = 0;
  lsv_cluster_max_multiplicity = 0;
  lsv_cluster_height = 0;
  
  fChain->SetBranchAddress("run", &run);
  fChain->SetBranchAddress("event_number", &event_number);
  fChain->SetBranchAddress("pps_counter", &pps_counter);
  fChain->SetBranchAddress("gps_fine_time_counter", &gps_fine_time_counter);
  fChain->SetBranchAddress("bad_time_alignment", &bad_time_alignment);
  fChain->SetBranchAddress("lsv_total_spe_charge", &lsv_total_spe_charge);
  fChain->SetBranchAddress("wt_total_spe_charge", &wt_total_spe_charge);
  fChain->SetBranchAddress("lsv_n_clusters", &lsv_n_clusters);
  fChain->SetBranchAddress("lsv_cluster_fixed_width_charge", &lsv_cluster_fixed_width_charge);
  fChain->SetBranchAddress("lsv_cluster_sat_corr_charge", &lsv_cluster_sat_corr_charge);
  fChain->SetBranchAddress("lsv_cluster_start_ns", &lsv_cluster_start_ns);
  fChain->SetBranchAddress("lsv_cluster_width_ns", &lsv_cluster_width_ns);
  fChain->SetBranchAddress("lsv_cluster_max_multiplicity", &lsv_cluster_max_multiplicity);
  fChain->SetBranchAddress("lsv_cluster_height", &lsv_cluster_height); 
}

int vetopulsesplit::Main_Fit(int start, int end, string Time)
{
  const int NUM=7;
  string Source[NUM] = {"C14","Co60","Co57","K40","Tl208","Th232","U238"};
  int NBins[NUM] = {200,3000,300,5000,2000,2000,3000};
  //  string MCinputdir = "/home/hqian/montecarlo/g4ds10/Linux-g++/pulse_splitter/out";
  string MCinputdir=Getinputdir();
  string last = "_cylinder_clustered.root";  
  //  vector<TH1F*> Hlist = {C14_EnergyMC,Co60_EnergyMC,Co57_EnergyMC,K40_EnergyMC,Tl208_EnergyMC,Th232_EnergyMC,U238_EnergyMC};  
  TObjArray Histlist(0);
  Histlist.Add(C14_EnergyMC);
  Histlist.Add(Co60_EnergyMC);
  Histlist.Add(Co57_EnergyMC);
  Histlist.Add(K40_EnergyMC);
  Histlist.Add(Tl208_EnergyMC);
  Histlist.Add(Th232_EnergyMC);
  Histlist.Add(U238_EnergyMC);
  for (int i=0; i<NUM; i++)
    {
      string MCinputfile = MCinputdir+"out" + Source[i].c_str() + last;
      TFile *file = new TFile(MCinputfile.c_str());
      TTree *dstree = (TTree *) file->Get("Recon");
      if(!file->IsOpen())
	{ cout<<"Error: Cannot Open MC Data File =>"<<MCinputfile<<endl;
	  return -1;
	}
      else  cout<<"Processing MC Data file: "<<MCinputfile<<endl;

      string temp = "_EnergyMC";
      string title = Source[i].c_str()+temp;
      //      Hlist.push_back(new TH1F(title.c_str(),title.c_str(),NBins,0,NBins*1.0));
      //      Hlist.at(i)->GetXaxis()->SetTitle("Energy [keV]");
      TH1F *temphist = dynamic_cast<TH1F*>(Histlist.At(i));
      temphist->SetTitle(title.c_str());
      temphist->SetBins(NBins[i],0,NBins[i]*1.0);
      Calculate_MCData(dstree,temphist);
    }
  
  Calculate_Norm(C14_EnergyMC,c14_energy_spectrum);
  Calculate_Convolution(C14_EnergyMC,C14_EnergyMC,C14_EnergyMC_1st);
  Calculate_Norm(C14_EnergyMC_1st,c14_energy_spectrum_1st);
  Calculate_Convolution(C14_EnergyMC,C14_EnergyMC_1st,C14_EnergyMC_2nd);
  Calculate_Norm(C14_EnergyMC_2nd,c14_energy_spectrum_2nd);
  Calculate_Convolution(C14_EnergyMC,C14_EnergyMC_2nd,C14_EnergyMC_3rd);
  Calculate_Norm(C14_EnergyMC_3rd,c14_energy_spectrum_3rd);

  Histlist.Add(C14_EnergyMC_1st);
  Histlist.Add(C14_EnergyMC_2nd);
  
  TChain *odtree =new TChain("odtree");
  Readdatafile(odtree,start,end);
  Int_t nEntries = odtree->GetEntries();
  cout<<"nEntries= "<<nEntries<<endl;

  double startrange=0;
  double endrange = 2100;
  TH1F* FullSpectrum = new TH1F("FullSpectrum","lsv cluster fixed width charge full spectrum;charge [PE];Counts",1000,startrange,endrange);
  TCut TimingCut = "lsv_cluster_start_ns.fArray>3770 && lsv_cluster_start_ns.fArray<3786";
  TCanvas *c1 = new TCanvas("c1", "Full Spectrum with Fit", 1200, 600);
  c1->SetLogy();
  c1->cd();
  odtree->Draw("lsv_cluster_fixed_width_charge.fArray>>FullSpectrum", TimingCut);
  FullSpectrum->Rebin(4);
  FullSpectrum->Sumw2();
  FullSpectrum->Scale(1.0/nEntries);

  const Int_t parnums =17;
  Double_t new_par[parnums];
  Double_t new_par_error[parnums];

  Double_t binw             = FullSpectrum->GetBinWidth(1);
  Double_t kB               = 0.00835;
  Double_t ly_mean          = 0.55;
  Double_t spe_var          = 0.14;
  Double_t ly_var           = 0.00642;
  Double_t constant         = 1.;
  Double_t basel_mean       = 0.386;
  Double_t basel_var        = 1.;
  Double_t threshold        = 1;
  Double_t window           = 2e-7;
  Double_t integralsum      = FullSpectrum->Integral();

  Double_t parvaluesamples[]={binw,kB,ly_mean,spe_var,ly_var,constant,basel_mean,basel_var,threshold,window,integralsum};
			      //,integralsum[1],integralsum[2],integralsum[3],integralsum[4],integralsum[5],integralsum[6]};

  string parnamesamples[]={"Bin Width","Birks' Constant","LY Mean [PE/keV]","Rel SPE Var","Rel LY Var",
			   "Constant","Baseline Mean [p.e]","Baseline Var","Threshold","Window Width [s]",
			   "C14 Amplitude","Co60 Amplitude","Co57 Amplitude","K40 Amplitude","Tl208 Amplitude","Th232 Amplitude","U238 Amplitude"};
  //  Double_t parlowlimits[parnums] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  //  Double_t paruplimits[]  = {1,10,0.01,0.7,0.14,1000,500,100,500,500,1e-6,1,1,1,1,1,1}
  Fit_C14 = new TF1("Fit_C14", C14Fit, startrange, endrange, parnums);           // Single C14
  Fit_C14_1st = new TF1("Fit_C14_1st",C14Fit_1st , startrange, endrange, parnums);   // 1st Pile-Up
  Fit_C14_2nd = new TF1("Fit_C14_2nd", C14Fit_2nd, startrange, endrange, parnums);   // 2nd Pile-Up
  //   FitC14 = new TF1("FitC14",C14Fit,start[0],end[0],parnums); //total C14 fit function
  //  FitC14->SetParNames(par0.c_str(),par1.c_str(),par2.c_str(),par3.c_str(),par4.c_str(),par5.c_str(),par6.c_str(),par7.c_str(),par8.c_str(),par9.c_str(),par10.c_str());
  
  Fit_Co60  = new TF1("Fit_Co60",Co60Fit,startrange,endrange,parnums);
  Fit_Co57  = new TF1("Fit_Co57",Co57Fit,startrange,endrange,parnums);
  Fit_K40   = new TF1("Fit_K40",K40Fit,startrange,endrange,parnums);
  Fit_Tl208 = new TF1("Fit_Tl208",Tl208Fit,startrange,endrange,parnums);
  Fit_Th232 = new TF1("Fit_Th232",Th232Fit,startrange,endrange,parnums);
  Fit_U238  = new TF1("Fit_U238",U238Fit,startrange,endrange,parnums);
  Fit_Total  = new TF1("Fit_Total",TotalFit,startrange,endrange,parnums);

  TObjArray Fitlist(0);
  Fitlist.Add(Fit_C14);
  Fitlist.Add(Fit_C14_1st);
  Fitlist.Add(Fit_C14_2nd);
  Fitlist.Add(Fit_Co60);
  Fitlist.Add(Fit_Co57);
  Fitlist.Add(Fit_K40);
  Fitlist.Add(Fit_Tl208);
  Fitlist.Add(Fit_Th232);
  Fitlist.Add(Fit_U238);
  Fitlist.Add(Fit_Total);
  
  vector<int> linecolor = Colors();
  for(int i=0; i<parnums; i++)
    {
      Fit_Total->SetParName(i,parnamesamples[i].c_str());
      if(i == 0 || i == 1 || i == 3) Fit_Total->FixParameter(i,parvaluesamples[i]);
      else{
	Fit_Total->SetParameter(i, parvaluesamples[i]);
	//      FitTotal->SetParLimits(i,parlowlimits[i],paruplimits[i]);
      }
    }
  Fit_Total->SetNpx(1000);
  //  FitTotal->SetLineColor(kRed);
  FullSpectrum->Fit(Fit_Total,"RV");
  Double_t chi2 = Fit_Total->GetChisquare();
  Fit_Total->GetParameters(new_par);

  for(Int_t i=0; i<Fitlist.GetSize(); i++)
    {
      TF1* tempfit = dynamic_cast<TF1*>(Fitlist.At(i));
      cout<<tempfit->GetName()<<endl;
      tempfit->SetLineColor(linecolor.at(i));
      tempfit->SetParameters(new_par);
      tempfit->Draw("SAME");
    }

  //  string outdir = "/home/hqian/montecarlo/g4ds10/Linux-g++/pulse_splitter/";
  string outdir=Getoutputdir();
  string output = outdir + "Pulse_Splitter_"+Time+".root";
  TFile outfile(output.c_str(), "RECREATE");
  Histlist.Write();
  Fitlist.Write();
  FullSpectrum->Write();
  c1->Write();
  outfile.Write();
  outfile.Close();

  std::cout << "==> Application finished." << std::endl;
  return 1;
}
  
vector<double> Energy_Spectrum_Data_row;
TH1F *temphist;
temphist = dynamic_cast<TH1F*>(Hlist.At(i));
cout<<temphist->GetTitle()<<endl;
Calculate_MCData(dstree,temphist, Energy_Spectrum_Data_row,i);
    Calculate_MCData(dstree,dynamic_cast<TH1F*>(Hlist.At(i)), Energy_Spectrum_Data_row);
  C14,Co60,Co57,K40,Tl208,Th232,U238
  setup the ranges for amplitude and ranges to fit with
Double_t lowrange[NUM] = {30,250,30,30,1150,0,880};
  Double_t uprange[NUM] = {80,1799.,80,810,1849,50,1300.};
  Int_t Histlowbin, Histupbin;
  Double_t integralsum[NUM];
  for (int i=0; i<NUM; i++)
    {
      Histlowbin = FullSpectrum->GetXaxis()->FindBin(lowrange[i]);
      Histupbin = FullSpectrum->GetXaxis()->FindBin(uprange[i]);
      integralsum[i] = FullSpectrum->Integral(Histlowbin,Histupbin);
      cout<<"integralsum for "<<Source[i]<<"= "<<integralsum[i]<<endl;
    }
  Double_t integralsample = FullSpectrum->Integral();
  Double_t start[NUM];
  Double_t end[NUM];
  for(int i=0; i<NUM; i++){
    start[i]=0;
    end[i]=1600;
  }
  //  FuncC14->SetLineColor(kPink);
  //  FuncPile1->SetLineColor(kOrange);
  //  FuncPile2->SetLineColor(kCyan);
  //  Int_t Colors[]={43,kMagenta,kBlue,kGreen,kYellow,29,7};
  
  TCanvas *c2 = new TCanvas("c2","C14 Convolution Energy Spectrum",1000,600);
  c2->Divide(2,2);
  c2->cd(1);
  C14_EnergyMC->Draw();
  c2->cd(2);
  C14_Energy_1st->Draw();
  c2->cd(3);
  C14_Energy_2nd->Draw();
  c2->cd(4);
  C14_Energy_3rd->Draw();
  TCanvas *c3 = new TCanvas("c3","C14 Fitter Function",1000,800);
  c3->Divide(2,2);
  c3->cd(1);
  FuncC14->Draw();
  c3->cd(2);
  FuncPile1->Draw();
  c3->cd(3);
  FuncPile2->Draw();
  c3->cd(4);
  FitC14->Draw();

  const int NUM=7;
  string Source[NUM] = {"C14","Co60","Co57","K40","Tl208","Th232","U238"};
  int NBins[NUM] = {200,3000,300,5000,2000,2000,3000};
  //  string MCinputdir = "/home/hqian/montecarlo/g4ds10/Linux-g++/pulse_splitter/out";
  string MCinputdir=Getinputdir();
  string last = "_cylinder_clustered.root";
  //  vector<TH1F*> Hlist = {C14_EnergyMC,Co60_EnergyMC,Co57_EnergyMC,K40_EnergyMC,Tl208_EnergyMC,Th232_EnergyMC,U238_EnergyMC};
  TObjArray Histlist(0);
  Histlist.Add(C14_EnergyMC);
  Histlist.Add(Co60_EnergyMC);
  Histlist.Add(Co57_EnergyMC);
  Histlist.Add(K40_EnergyMC);
  Histlist.Add(Tl208_EnergyMC);
  Histlist.Add(Th232_EnergyMC);
  Histlist.Add(U238_EnergyMC);
  for (int i=0; i<NUM; i++)
  {
  string MCinputfile = MCinputdir+"out" + Source[i].c_str() + last;
  TFile *file = new TFile(MCinputfile.c_str());
  TTree *dstree = (TTree *) file->Get("Recon");
  if(!file->IsOpen())
  { cout<<"Error: Cannot Open MC Data File =>"<<MCinputfile<<endl;
  mp;
  
  split->Calculate_Convolution(C14_EnergyMC,C14_EnergyMC_1st,C14_Ene
  Histlist.Add(C14_EnergyMC_1st);                                                                                            
  
  Histlist(0);
  Histlist->Add(C14_EnergyMC);
  cout<<"2"<<endl;
  Histlist->Add(Co60_EnergyMC);
  Histlist->Add(Co57_EnergyMC);
  Histlist->Add(K40_EnergyMC);
  Histlist->Add(Tl208_EnergyMC);
  Histlist->Add(Th232_EnergyMC);
  Histlist->Add(U238_EnergyMC);

  Calculate_Norm(C14_EnergyMC,c14_energy_spectrum);
  Calculate_Convolution(C14_EnergyMC,C14_EnergyMC,C14_EnergyMC_1st);
  Calculate_Norm(C14_EnergyMC_1st,c14_energy_spectrum_1st);
  Calculate_Convolution(C14_EnergyMC,C14_EnergyMC_1st,C14_EnergyMC_2nd);
  Calculate_Norm(C14_EnergyMC_2nd,c14_energy_spectrum_2nd);
  Calculate_Convolution(C14_EnergyMC,C14_EnergyMC_2nd,C14_EnergyMC_3rd);
  Calculate_Norm(C14_EnergyMC_3rd,c14_energy_spectrum_3rd);

  //      string title = Source[i].c_str()+temp;
  //     TH1F *temphist = dynamic_cast<TH1F*>(Histlist->At(i));
  //      temphist->SetTitle(title.c_str());
  //      temphist->SetBins(NBins[i],0,NBins[i]*1.0);
  Histlist->Add(C14_EnergyMC_1st);
  Histlist->Add(C14_EnergyMC_2nd);
  Histlist->Add(C14_EnergyMC_3rd);
  
  }
  else{
    for(Int_t i=0; i<NEntries; i++)
      {
	sum = 0;
	dstree->GetEntry(i);
	vector<double> tagene, timeinlsv, timeintpc;
	for(Int_t j=0; j<et->size(); j++)
	  {
	    if(volume->at(j)=="p_active")
	      timeintpc.push_back(et->at(j)*1.e+9);
	    if(volume->at(j)=="p_scint")
	      {
		timeinlsv.push_back(et->at(j)*1.e+9);
		tagene.push_back(eqch->at(j));
	      }
	  }
	if(!timeintpc.empty() && !timeinlsv.empty())
	  {
	    for(int k=0; k<timeinlsv.size(); k++){
	      double time_diff = timeinlsv.at(k) - timeintpc.at(0);
	      if(time_diff<8. && time_diff>0)
		sum += tagene.at(k);
	    }
	  }
	//      Energy_Spectrum.push_back(sum);
	Energy_Spectrum->Fill(sum);
      }
  }

  for(Int_t i=1; i<=NBins; i++)
    {
      //      Energy_new->Fill(C14_Energy->GetBinCenter(i),C14_Energy->GetBinContent(i));
      Energy_Spectrum_Data.push_back(Energy_Spectrum->GetBinContent(i));
    }
  
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

double vetopulsesplit::C14Fit(double* x, double* params) {
  int np = c14_energy_spectrum.size();
  double q = x[0];
  double de = (c14_spectrum_endpoint - c14_spectrum_startpoint)/(np);
  double result = 0;
  int choice = 1;

  for (int i = 0; i < np; i++) {
    double energy = de*i;
    int factor = ((i % 2) ? 4:2);
    if (i == 0 || i + 1 == np) factor = 1;
    result += factor*(c14_energy_spectrum.at(i)*response_function(q, energy, params, choice));
  } result *= de/3;

  double C14Rate = params[9];
  return result*(C14Rate)*params[0];
  //  return result*(C14Rate)*params[0]*BGFraction[0];

}

//// 1ST PILE-UP FUNCTION  ////////////////////////////////////////////////////////////
double vetopulsesplit::C14Fit_1st (double* x, double* params) {
  int np = c14_energy_spectrum_1st.size();
  double q = x[0];
  double de = (2*c14_spectrum_endpoint- c14_spectrum_startpoint)/(np);
  double result = 0;
  int choice = 0; // 1st pile-up
  double window = params[8];

  for (int i = 0; i < np; i++) {
    double energy = de*i;
    int factor = ((i % 2) ? 4:2);
    if (i == 0 || i + 1 == np) factor = 1;
    result += factor*(c14_energy_spectrum_1st.at(i)*response_function(q, energy, params, choice));
  } result *= de/3;

  //double stRate = params[8]*params[8]*window;
  double stRate = pow(params[9],2)*window*exp(-params[9]*window);
  return result*(stRate)*params[0];
  //return result*(stRate)*params[0]*BGFraction[0];

}

//// 2ND PILE-UP FUNCTION  ////////////////////////////////////////////////////////////
double vetopulsesplit::C14Fit_2nd (double* x, double* params) {
  int np = c14_energy_spectrum_2nd.size();
  double q = x[0];
  double de = (3*c14_spectrum_endpoint- c14_spectrum_startpoint)/(np);
  double result = 0;
  int choice = 0; // 2nd pile-up
  double window = params[8];

  for (int i = 0; i < np; i++) {
    double energy = de*i;
    int factor = ((i % 2) ? 4:2);
    if (i == 0 || i + 1 == np) factor = 1;
    result += factor*(c14_energy_spectrum_2nd.at(i)*response_function(q, energy, params, choice));
  } result *= de/3;

  //double ndRate = params[8]*params[8]*params[8]*window*window;
  double ndRate = pow(params[9],3)*pow(window,2)*0.5*exp(-params[9]*window);
  return result*(ndRate)*params[0];
  //  return result*(ndRate)*params[0]*BGFraction[0];

}

 
  const int NUM=7;
  const int piles=2;
  string Source[NUM] = {"C14","Co60","Co57","K40","Tl208","Th232","U238"};
  int NBins[NUM] = {160,2600,800,1600,4500,3500,3500};
  string Source[NUM] = {"C14","Co60","Co57","K40","Th232","U238","U235"};
  int NBins[NUM] = {160,2600,800,1600,3500,3500,1000};
  EnergyMC = new TH1F*[NUM];
  VetoMC_ntuple = new TNtuple*[NUM];
 

void vetopulsesplit::Recon_DataAnalysis()
{
  SetIsotopes();
  string MCinputdir=GetMCinputdir();
  string label = "_v1";
  //  string label = "";
  string last = label+"_clustered.root";  
  string temp = "_EnergyMC";
  for(size_t i=0; i<Source.size(); i++)  
    //for(size_t i=0; i<2; i++)  
    {
  
      string MCinputfile = MCinputdir+"out" + Source[i]+ last;
      TFile *file = new TFile(MCinputfile.c_str());
      TTree *dstree = (TTree *) file->Get("Recon");
      if(!file->IsOpen())
        cout<<"Error: Cannot Open MC Data File =>"<<MCinputfile<<endl;       
      else  cout<<"Processing MC Data file: "<<MCinputfile<<endl;
  
      TString hname=Form("%s%s",Source[i].c_str(),temp.c_str());
      TString htitle=Form("%s MC Energy",Source[i].c_str());
      EnergyMC.push_back(new TH1F(hname,htitle,NBins[i],0,NBins[i]*1.0));
      Histlist.Add(EnergyMC[i]);      
  
      //      string tag = "out"+Source.at(i)+label+".root";
      Dstree_Analysis(label,i);
      if(i==0)
	//	fraction.push_back(Calculate_MCData(dstree,EnergyMC[i]));
	fraction.push_back(Calculate_MCData(i,last,EnergyMC[i]));
      else
	{
	TString hname_ntuple=Form("%s_vetoMCntuple",Source[i].c_str());             
	  VetoMC_ntuple.push_back(new TNtuple(hname_ntuple,Source[i].c_str(),"eqch:edep:gps:et"));
	  Histlist.Add(VetoMC_ntuple.back());	  
	  Calculate_MCCoincidenceData(dstree,EnergyMC[i],VetoMC_ntuple.back());
	
	  //  fraction.push_back(Calculate_MCCoincidenceData(dstree,EnergyMC[i]));    
	  fraction.push_back(Calculate_MCCoincidenceData(i,last,EnergyMC[i]));    
	}
      cout<<fraction.back()<<endl;
      //      dstree->SetDirectory(0);      
    }  

  C14_EnergyMC_1st = new TH1F("C14_EnergyMC_1st","1st Convolution of C14 MC Energy",NBins[0]*2,0,NBins[0]*2.0);
  C14_EnergyMC_2nd = new TH1F("C14_EnergyMC_2nd","2nd Convolution of C14 MC Energy",NBins[0]*3,0,NBins[0]*3.0);
  C14_EnergyMC_3rd = new TH1F("C14_EnergyMC_3rd","3rd Convolution of C14 MC Energy",NBins[0]*4,0,NBins[0]*4.0);
  Calculate_Convolution(EnergyMC[0],EnergyMC[0],C14_EnergyMC_1st);
  Calculate_Convolution(EnergyMC[0],C14_EnergyMC_1st,C14_EnergyMC_2nd);
  Calculate_Convolution(EnergyMC[0],C14_EnergyMC_2nd,C14_EnergyMC_3rd);  
  Histlist.Add(C14_EnergyMC_1st);
  Histlist.Add(C14_EnergyMC_2nd);
  Histlist.Add(C14_EnergyMC_3rd);

  string outdir=GetMCinputdir();
  string output = outdir +"PulseSplitMCEnergy"+label+".root";
  //  string output = outdir +"PulseSplitMCEnergy.root";
  TFile outfile(output.c_str(), "RECREATE");
  Histlist.Write();                       
  outfile.Write();
  outfile.Close();
  cout<<"Successfully save MC Data."<<endl;
}  

string MCinputdir=GetMCinputdir();  
string MCinputfile = MCinputdir+"out" + Source.at(k)+ last;
TFile *file = new TFile(MCinputfile.c_str());
TTree *dstree = (TTree *) file->Get("Recon");
if(!file->IsOpen())
cout<<"Error: Cannot Open MC Data File =>"<<MCinputfile<<endl;       
else  cout<<"Processing MC Data file: "<<MCinputfile<<endl;  


string MCinputdir=GetMCinputdir();
  string MCinputfile = MCinputdir+"out"+Source.at(k)+label+".root";
  TFile *file = new TFile(MCinputfile.c_str());
  TTree *dstree = (TTree *) file->Get("dstree");
  if(!file->IsOpen())
    cout<<"Error: Cannot Open MC Data File =>"<<MCinputfile<<endl;
  else  cout<<"Processing MC Data file: "<<MCinputfile<<endl;


  Int_t           ev;
  Float_t         tpcene;
  Float_t         vetoene;
  Float_t         ene;
  Int_t           ndeposits;
  Int_t           dep_pdg[1103];   //[ndeposits]
  Int_t           dep_mat[1103];   //[ndeposits]
  Int_t           dep_track[1103];   //[ndeposits]
  Int_t           dep_parenttrack[1103];   //[ndeposits]
  Double_t        dep_time[1103];   //[ndeposits]
  Float_t         dep_ene[1103];   //[ndeposits]
  Float_t         dep_x[1103];   //[ndeposits]
  Float_t         dep_y[1103];   //[ndeposits]
  Float_t         dep_z[1103];   //[ndeposits]      

  dstree->SetBranchAddress("ev", &ev);
  dstree->SetBranchAddress("tpcene", &tpcene);
  dstree->SetBranchAddress("vetoene", &vetoene);
  dstree->SetBranchAddress("ene", &ene);
  dstree->SetBranchAddress("ndeposits", &ndeposits);
  dstree->SetBranchAddress("dep_pdg", dep_pdg);
  dstree->SetBranchAddress("dep_mat", dep_mat);
  dstree->SetBranchAddress("dep_track", dep_track);
  dstree->SetBranchAddress("dep_parenttrack", dep_parenttrack);
  dstree->SetBranchAddress("dep_time", dep_time);
  dstree->SetBranchAddress("dep_ene", dep_ene);
  dstree->SetBranchAddress("dep_x", dep_x);
  dstree->SetBranchAddress("dep_y", dep_y);
  dstree->SetBranchAddress("dep_z", dep_z);


double vetopulsesplit::Ni65Fit (double* x, double* params){
  int nBins = EnergyMC[8]->GetNbinsX();
  double delta = EnergyMC[8]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[8]->GetBinCenter(i);
      spectrum = EnergyMC[8]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }

  double Ni65Rate = params[17];
  return result*(Ni65Rate)*params[0]*delta;

}

double vetopulsesplit::Fe59Fit (double* x, double* params){
  int nBins = EnergyMC[9]->GetNbinsX();
  double delta = EnergyMC[9]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[9]->GetBinCenter(i);
      spectrum = EnergyMC[9]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }

  double Fe59Rate = params[18];
  return result*(Fe59Rate)*params[0]*delta;

}


void vetopulsesplit::SetIsotopes()
{
  Source.push_back("C14");
  Source.push_back("Po210");
  Source.push_back("Bi210");
  Source.push_back("Rn212");
  Source.push_back("Po212");
  Source.push_back("Bi212");
  Source.push_back("Rn214");
  Source.push_back("Po214");
  Source.push_back("Bi214");
  Source.push_back("Bi214stem");
  Source.push_back("Tl208");
  Source.push_back("K40");
  Source.push_back("Co60");
  
  PDG.push_back(1000060140);
  PDG.push_back(1000842100);
  PDG.push_back(1000832100);
  PDG.push_back(1000862120);
  PDG.push_back(1000842120);
  PDG.push_back(1000832120);
  PDG.push_back(1000862140);
  PDG.push_back(1000842140);
  PDG.push_back(1000832140);
  PDG.push_back(1000832140);
  PDG.push_back(1000812080);
  PDG.push_back(1000270600);
  
  NBins.push_back(160);
  NBins.push_back(160);
  NBins.push_back(160);
  NBins.push_back(160);
  NBins.push_back(160);
  NBins.push_back(160);
  NBins.push_back(160);
  NBins.push_back(160);
  NBins.push_back(160);
  NBins.push_back(160);
 
}
*/






