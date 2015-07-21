#include "c14pulse.hh"
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>

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

using namespace std;
ClassImp(c14pulse);

double c14pulse::Isotope_Sum(int isotopes, double *isotope_activity)
{
  double activity_sum=0;
  for(int i=0; i<isotopes; i++)    
    activity_sum += isotope_activity[i];
  return activity_sum;
}

void c14pulse::Isotope_Fraction(int isotopes, double *isotope_activity, double *isotope_fraction)
{
  double sum = Isotope_Sum(isotopes,isotope_activity);
  for(int i=0; i<isotopes; i++)   
    {
      isotope_fraction[i] = isotope_activity[i]*1.0/sum;
      BGFraction.push_back(isotope_fraction[i]);
    }
}

void c14pulse::Calculate_MCData(TTree *dstree, TH1F *Energy_Spectrum)
{
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
  
  for(Int_t i=0; i<NEntries; i++)
    { dstree->GetEntry(i);
      for(Int_t j=0; j<et->size(); j++)
	{
	  if(volume->at(j)=="p_scint")
	    Energy_Spectrum->Fill(eqch->at(j));
	}
    }
  Int_t NBins = Energy_Spectrum->GetNbinsX();
  cout<<"NBins= "<<NBins<<endl;
  Double_t integralsum = Energy_Spectrum->Integral();
  Energy_Spectrum->Scale(1./integralsum);
}

bool c14pulse::multicut(float height,float multiplicity, float charge){
  return height/multiplicity < (2.563e7 + TMath::Sqrt(1.574e14+1.390e12*(charge-14.40)*(charge-14.40)));
}

void c14pulse::Readdatafile(TChain *t,int i)
{
  //  string dirname="/darkside/users/hqian/pulse_splitter/";
  string dirname=GetRealinputdir();
  string filename;
  string middle="ODRun000";
  string last=".root";
  stringstream oss;
  ifstream NameCheck;

  oss<<i;
  filename=dirname+middle+oss.str()+last;
  NameCheck.open(filename.c_str());
  if(!NameCheck.good())
    {
      oss.str("");
      NameCheck.close();
      cout<<"ERROR: File does not exist!"<<endl;
      //      return false;
    }
  else{
    t->Add(filename.c_str());
    cout<<"Processing Data file: "<<filename<<endl;
    oss.str("");
  }
  NameCheck.close();
  //  return true;
}

void c14pulse::Readdatafile(TChain *t, int startfile, int endfile)
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

void c14pulse::Calculate_Convolution(TH1F *energy_spectrum1, TH1F *energy_spectrum2, TH1F *convolution)
{
  int np1 = energy_spectrum1->GetNbinsX();
  int np2 = energy_spectrum2->GetNbinsX();
  int nconv = np1 + np2;
  double binsize;
  cout<<"Calculate Convolution!!"<<endl;
  for (int n = 1; n <= nconv; n++) {
    double sum = 0;
    for (int m = 1; m <= n; m++) {
      if (m > np1  || (n - m) > np2 ) continue;
      sum += (energy_spectrum1->GetBinContent(m)) * (energy_spectrum2->GetBinContent(n-m+1));
    }
    convolution->SetBinContent(n,sum);
   }
  cout<<"Finish Calculating Convolution!!"<<endl;
}

void c14pulse::Calculate_Norm (TH1F *energy_spectrum,  vector<double> &convolution_spectrum)
{
  int np = energy_spectrum->GetNbinsX();
  double integral = energy_spectrum->Integral();
  cout<<"integral= "<<integral<<endl;
  energy_spectrum->Scale(1./integral);
  for(int i=1; i<=np; i++)
    convolution_spectrum.push_back(energy_spectrum->GetBinContent(i));
}

double c14pulse::scint_e_quenching (double e_keV, double* params) {

  double kB =0.012;// params[1];
  double A1 = -0.6292-0.2181*log(kB);
  double A2 = -0.3057-0.1024*log(kB);
  double A3 = -0.0673-0.03353*log(kB);
  double A4 = 0.009876+0.002276*log(kB);
  double A5 = -0.2814-0.09964*log(kB);
  double A6 = -0.09962-0.0376*log(kB);

  
  return ( (A1 + A2*log(e_keV) + A3*log(e_keV)*log(e_keV) + A4*log(e_keV)*log(e_keV)*log(e_keV))/(1 + A5*log(e_keV) + A6*log(e_keV)*log(e_keV) + A4*log(e_keV)*log(e_keV)*log(e_keV)) );
}

double c14pulse::scint_e_quenching (double e_keV) {

  double kB =0.012;// params[1];
  double A1 = -0.6292-0.2181*log(kB);
  double A2 = -0.3057-0.1024*log(kB);
  double A3 = -0.0673-0.03353*log(kB);
  double A4 = 0.009876+0.002276*log(kB);
  double A5 = -0.2814-0.09964*log(kB);
  double A6 = -0.09962-0.0376*log(kB);
  /*
  double A1 = 0.32903;
  double A2 = 0.14404;
  double A3 = 0.08059;
  double A4 = -0.04536E-3;
  double A5 = 0.15623;
  double A6 = 0.06611;
  */
  return ( (A1 + A2*log(e_keV) + A3*log(e_keV)*log(e_keV) + A4*log(e_keV)*log(e_keV)*log(e_keV))/(1 + A5*log(e_keV) + A6*log(e_keV)*log(e_keV) + A4*log(e_keV)*log(e_keV)*log(e_keV)) );
}

void c14pulse::Scint_Quenching_Model(TTree *dstree, TH1F *scint_quench_model, TH1F *quenching_plot)
{
  Int_t NEntries = dstree->GetEntries();
  cout<<"NEntries= "<<NEntries<<endl;

  /*  vector<TString> *volume;
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
  
  for(Int_t i=0; i<NEntries; i++)
    { dstree->GetEntry(i);
      for(Int_t j=0; j<et->size(); j++)
	{
	  if(volume->at(j)=="p_scint")	    
	    {
	      double quench=scint_e_quenching(edep->at(j));
	      //     scint_quench_model->Fill(edep->at(j),scint_e_quenching(edep->at(j)));
	      int bin = scint_quench_model->FindBin(edep->at(j));
	      scint_quench_model->SetBinContent(bin,quench);
	      quenching_plot->Fill(edep->at(j)*(scint_e_quenching(edep->at(j))));
	    }
	}
    }
  Int_t NBins = quenching_plot->GetNbinsX();
  cout<<"NBins= "<<NBins<<endl;
*/
  const Int_t N = 100000;
  Float_t   dep_ene[N];
  Int_t     nDeposits;
  Int_t     dep_mat[N];
  dstree->SetBranchAddress("dep_ene",dep_ene);
  dstree->SetBranchAddress("ndeposits",&nDeposits);
  dstree->SetBranchAddress("dep_mat",dep_mat);

  for(Int_t i=0; i<NEntries; i++)
    {
      dstree->GetEntry(i);
      double edep_sum=0;
      double eqch_sum=0;
      for(Int_t j=0; j<nDeposits; j++)
	{
	  if(dep_mat[j] == 2)
	    {
	      double edep = dep_ene[j]*1000.;
	      double quench=scint_e_quenching(edep);
	      double eqch= edep*quench;
	      edep_sum += edep;
	      eqch_sum += eqch;
	    }}
      quenching_plot->Fill(eqch_sum);
      scint_quench_model->Fill(edep_sum);      
    }
  quenching_plot->Scale(1./(quenching_plot->Integral()));
  scint_quench_model->Scale(1./(scint_quench_model->Integral()));
}


double c14pulse::response_function (double q, double energy, double* params, int choice) {

  double ly_mean = params[1];
  double spe_var = params[2];
  double ly_var = params[3];
  double baseline_mean = params[5];
  double baseline_var = params[6];
  //  double threshold = params[8];
  double pe;
  double gaus_mean;

  pe = energy*ly_mean;
  //  if (energy == 0) pe = energy*ly_mean;
  //pe = energy*ly_mean*scint_e_quenching(energy, params);

  if (choice == 1) gaus_mean = pe + (baseline_mean + params[7]);
  else gaus_mean = pe + baseline_mean;

  double gaus_var = baseline_var + (1 + spe_var)*pe + ly_var*pe*pe;
  //  double gaus_var = baseline_var*(1+spe_var)*gaus_mean+ly_var*TMath::Power(gaus_mean,2);
  double gaus_var_inv = 1.0 / gaus_var;
  double arg = -0.5*(q - gaus_mean)*(q - gaus_mean)*gaus_var_inv;

  return 0.3989422804 * sqrt(gaus_var_inv) * exp(arg);
}
/*
double c14pulse::C14Fit(double* x, double* params) {
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
  //  return result*(C14Rate)*params[0];
  return result*(C14Rate)*params[0];

}

//// 1ST PILE-UP FUNCTION  ////////////////////////////////////////////////////////////
double c14pulse::C14Fit_1st (double* x, double* params) {
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
  //  return result*(stRate)*params[0];
  return result*(stRate)*params[0];

}

//// 2ND PILE-UP FUNCTION  ////////////////////////////////////////////////////////////
double c14pulse::C14Fit_2nd (double* x, double* params) {
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
  //  return result*(ndRate)*params[0];
  return result*(ndRate)*params[0];
}
*/
double c14pulse::C14Fit (double* x, double* params){
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
  return result*(C14Rate)*delta*params[0];
  //  return result*(C14Rate)*params[0];
  //  return result*(C14Rate)*params[0]*BGFraction[1];
}

double c14pulse::C14Fit_1st (double* x, double* params){
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
  return result*(stRate)*delta*params[0];
  //  return result*(stRate)*params[0];
  //  return result*(stRate)*params[0]*BGFraction[1];
}

double c14pulse::C14Fit_2nd (double* x, double* params){
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
  return result*(ndRate)*delta*params[0];
  //  return result*(ndRate)*params[0];
  //  return result*(ndRate)*params[0]*BGFraction[1];
}



/*
double c14pulse::Co60Fit (double* x, double* params){
  int nBins = EnergyMC[1]->GetNbinsX();
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
  //  return result*(Co60Rate)*params[0];
  return result*(Co60Rate)*params[0]*BGFraction[1];

}

double c14pulse::Co57Fit (double* x, double* params){
  int nBins = EnergyMC[2]->GetNbinsX();
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
  //  return result*(Co57Rate)*params[0];
  return result*(Co57Rate)*params[0]*BGFraction[2];

}

double c14pulse::K40Fit (double* x, double* params){
  int nBins = EnergyMC[3]->GetNbinsX();
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
  //  return result*(K40Rate)*params[0];
  return result*(K40Rate)*params[0]*BGFraction[3];

}

double c14pulse::Tl208Fit (double* x, double* params){
  int nBins = EnergyMC[4]->GetNbinsX();
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
  return result*(Tl208Rate)*params[0];
}

double c14pulse::Th232Fit (double* x, double* params){
  int nBins = EnergyMC[4]->GetNbinsX();
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
  //  return result*(Th232Rate)*params[0];
  return result*(Th232Rate)*params[0]*BGFraction[4];

}

double c14pulse::U238Fit (double* x, double* params){
  int nBins = EnergyMC[5]->GetNbinsX();
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
  double U238Rate = params[14];
  //  return result*(U238Rate)*params[0];
  return result*(U238Rate)*params[0]*BGFraction[5];

}
*/

double c14pulse::TotalFit(Double_t *x, Double_t *params)
{
  //  double C14_value = C14(x, params);
  // double stPile_value = stPile(x, params);
  // double ndPile_value = ndPile(x, params);
  double C14_value = C14Fit(x,params);
  double C14_1st_value = C14Fit_1st(x,params);
  double C14_2nd_value = C14Fit_2nd(x,params);
  /*  double Co60_value = Co60Fit(x,params);
  double Co57_value = Co57Fit(x,params);
  double K40_value = K40Fit(x,params);
  //  double Tl208_value = Tl208Fit(x,params);
  double Th232_value = Th232Fit(x,params);
  double U238_value = U238Fit(x,params);
  */
  //  return C14_value + stPile_value + ndPile_value + Co60_value+Co57_value+K40_value+Tl208_value+Th232_value+U238_value+ params[6];
  //  return C14_value +C14_1st_value+C14_2nd_value+ Co60_value+Co57_value+K40_value+Tl208_value+Th232_value+U238_value+params[4];
  //  return C14_value +C14_1st_value+C14_2nd_value+ Co60_value+Co57_value+K40_value+Th232_value+U238_value+params[4];
  /*  return C14_value*BGFraction[0] +C14_1st_value*BGFraction[0]+C14_2nd_value*BGFraction[0]
    + Co60_value*BGFraction[1]+Co57_value*BGFraction[2]+K40_value*BGFraction[3]
    +Th232_value*BGFraction[4]+U238_value*BGFraction[5]+params[4];
  */
  return C14_value + C14_1st_value +C14_2nd_value +params[4];
}

void c14pulse::Recon_DataAnalysis()
{
  const int NUM=6;
  //  const int piles=2;
  //  string Source[NUM] = {"C14","Co60","Co57","K40","Tl208","Th232","U238"};
  //int NBins[NUM] = {160,2600,800,1600,4500,3500,3500};
  string Source[NUM] = {"C14","Co60","Co57","K40","Th232","U238"};
  int NBins[NUM] = {160,2600,800,1600,3500,3500};
  string MCinputdir=GetMCinputdir();
  string last = "_v4_clustered.root";  
  //  string last = "_clustered.root";  
  EnergyMC = new TH1F*[NUM];
  C14_Quenching_Model = new TH1F("C14_Quenching_Model","C14 Quenching Model;edep [keV];quenching factor",NBins[0],0,NBins[0]);
  quenching_plot = new TH1F("quenching_plot","C14 Quenching Model Plot;eqch [keVee]",NBins[0],0,NBins[0]);
  for(int i=0; i<NUM-5; i++)
    {
      cout<<"3"<<endl;
      string MCinputfile = MCinputdir+"out" + Source[i]+ last;
      TFile *file = new TFile(MCinputfile.c_str());
      TTree *dstree = (TTree *) file->Get("Recon");
      if(!file->IsOpen())
        cout<<"Error: Cannot Open MC Data File =>"<<MCinputfile<<endl;       
      else  cout<<"Processing MC Data file: "<<MCinputfile<<endl;
      string temp = "_EnergyMC";
      TString hname=Form("%s%s",Source[i].c_str(),temp.c_str());
      TString htitle=Form("%s MC Energy",Source[i].c_str());
      EnergyMC[i] = new TH1F(hname,htitle,NBins[i],0,NBins[i]*1.0);
      Histlist.Add(EnergyMC[i]);
      Calculate_MCData(dstree,EnergyMC[i]);
      //  EnergyMC[i]->Sumw2();
    }  
  C14_EnergyMC_1st = new TH1F("C14_EnergyMC_1st","1st Convolution of C14 MC Energy",NBins[0]*2,0,NBins[0]*2.0);
  C14_EnergyMC_2nd = new TH1F("C14_EnergyMC_2nd","2nd Convolution of C14 MC Energy",NBins[0]*3,0,NBins[0]*3.0);
  C14_EnergyMC_3rd = new TH1F("C14_EnergyMC_3rd","3rd Convolution of C14 MC Energy",NBins[0]*4,0,NBins[0]*4.0);

  Calculate_Norm(EnergyMC[0],c14_energy_spectrum);
  Calculate_Convolution(EnergyMC[0],EnergyMC[0],C14_EnergyMC_1st);
  Calculate_Norm(C14_EnergyMC_1st,c14_energy_spectrum_1st);
  Calculate_Convolution(EnergyMC[0],C14_EnergyMC_1st,C14_EnergyMC_2nd);
  Calculate_Norm(C14_EnergyMC_2nd,c14_energy_spectrum_2nd);
  Calculate_Convolution(EnergyMC[0],C14_EnergyMC_2nd,C14_EnergyMC_3rd);
  Calculate_Norm(C14_EnergyMC_3rd,c14_energy_spectrum_3rd);
  
  Histlist.Add(C14_EnergyMC_1st);
  Histlist.Add(C14_EnergyMC_2nd);
  Histlist.Add(C14_EnergyMC_3rd);

  string MCinputfile = MCinputdir+"out" + Source[0]+ ".root";
  TFile *file = new TFile(MCinputfile.c_str());
  TTree *dstree = (TTree *) file->Get("dstree");
  Scint_Quenching_Model(dstree,C14_Quenching_Model,quenching_plot);
  Histlist.Add(C14_Quenching_Model);
  Histlist.Add(quenching_plot);
  
  /*  C14_EnergyMC_1st->Sumw2();
  C14_EnergyMC_2nd->Sumw2();
  C14_EnergyMC_3rd->Sumw2();
  */
  
  string outdir=GetMCinputdir();
  string output = outdir +"C14PulseSplitMCEnergy_v4.root";
  //  string output = outdir +"PulseSplitMCEnergy.root";
  TFile outfile(output.c_str(), "RECREATE");
  Histlist.Write();                       
  outfile.Write();
  outfile.Close();
  cout<<"Successfully save MC Data."<<endl;
}  

void c14pulse::ODTree_DataAnalysis(int start, int end)
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

  TH1F* FullSpectrum = new TH1F("FullSpectrum","lsv cluster fixed width charge full spectrum;charge [PE];Counts",Bins,0,2000);
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
void c14pulse::ODTree_Branch(TChain *fChain)
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
*/







/*
int c14pulse::Main_Fit(int start, int end, string Time)
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
*/
  
/*      vector<double> Energy_Spectrum_Data_row;
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
  */

 /*  const int NUM=7;
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
       
 */
    /* }
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
    */
  /*  for(Int_t i=1; i<=NBins; i++)
    {
      //      Energy_new->Fill(C14_Energy->GetBinCenter(i),C14_Energy->GetBinContent(i));
      Energy_Spectrum_Data.push_back(Energy_Spectrum->GetBinContent(i));
    }
  */
