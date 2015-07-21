#ifndef _vetopulsesplit_H
#define _vetopulsesplit_H
#include <map>

#include "TH1F.h"
#include "TH2F.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TRint.h"
#include "TColor.h"
#include "TNtuple.h"

using namespace std;
class vetopulsesplit:public TObject 
{
public :
  vetopulsesplit():
    C14_EnergyMC_1st(0),C14_EnergyMC_2nd(0), C14_EnergyMC_3rd(0),
    MCValue_ntuple(0)
    //    fraction(0)
  {// c14_spectrum_endpoint = 156.27; 
    //    c14_spectrum_startpoint = 0;     
  }
  virtual ~vetopulsesplit() { 
    MCFile.clear();
    MCinputdir.clear();
    Realinputdir.clear();
    outputdir.clear();
    C14_EnergyMC_1st=NULL;
    C14_EnergyMC_2nd=NULL;
    C14_EnergyMC_3rd=NULL;
  }
  void Recon_DataAnalysis();
  void ODTree_DataAnalysis(int,int);
  void Dstree_Analysis(int);

  //  double Calculate_MCData(TTree*,TH1F*);
  //  void Calculate_MCCoincidenceData(TTree*,TH1F*,TNtuple*);
  // double Calculate_MCCoincidenceData(TTree*,TH1F*);
  void Calculate_MCData(int,TH1F*);
  void Calculate_MCCoincidenceData(int,TH1F*);
  void Calculate_U238MCCoincidenceData(int,TH1F*,TH1F*);
  
  //  void Readdatafile(TChain*,int);
  void Readdatafile(TChain*, int, int);
  void DST_Readdatafile(TChain*, int, int);
  void Recon_Readdatafile(TChain*, int,int,int,string);

  void Calculate_Convolution(TH1F*, TH1F*, TH1F*);
  void Calculate_Norm(TH1F*, vector<double> &);
  double scint_e_quenching(double, double*);
  double response_function(double, double, double*,int);
  // int Main_Fit(int, int,string);
  double C14Fit(double*,double*);
  double C14Fit_1st(double*,double*);
  double C14Fit_2nd(double*,double*);
  double Co60Fit(double*,double*);
  double Co57Fit(double*,double*);
  double K40Fit(double*,double*);
  double Tl208Fit(double*,double*);
  double Th232Fit(double*,double*);
  double Th232LowerFit(double*,double*);
  double U238UpperFit(double*,double*);  
  double U238LowerFit(double*,double*);  
  double U238Fit(double*,double*);  
  double Rn222Fit(double*,double*);  
  double K42Fit(double*,double*);  
  double U235Fit(double*,double*);  
  double U235LowerFit(double*,double*);  
  double TotalFit(double*,double*);  
  TObjArray GetHistArray() { return Histlist; }
  string GetMCinputdir() { return MCinputdir; }
  void SetMCFile(string val) { MCFile=val; }
  string GetMCFile() { return MCFile; }
  void SetMCinputdir(string val) { MCinputdir=val; }
  string GetRealinputdir() { return Realinputdir; }
  void SetRealinputdir(string val) { Realinputdir=val; }
  string Getoutputdir() { return outputdir; }
  void Setoutputdir(string val) { outputdir=val; }
  double Isotope_Sum(int,vector<double>);
  void Isotope_Fraction(int,vector<double>,double*);
  bool multicut(float,float,float);
  void SetIsotopes();
  vector<string> GetIsotopes() { return Source; }
  vector<double> GetFraction() { return fraction; } 
  vector<int> GetNBins() { return NBins; }
  vector<int> Colors();
  vector<TH1F*>* GetEnergyMC() { return &EnergyMC; }
  void SetEnergyMC(vector<TH1F*> val) { EnergyMC=val; }
  

private:
  string MCFile;
  string MCinputdir;
  string Realinputdir;
  string outputdir;
  vector<string> Source;
  vector<Int_t> PDG;
  vector<int> NBins;
  vector<int> TBins;
  vector<int> Normalization;
  TObjArray Histlist;
  TH1F *C14_EnergyMC_1st, *C14_EnergyMC_2nd, *C14_EnergyMC_3rd;
  vector<TH1F*> EnergyMC;
  vector<TH1F*> EdepMC;
  vector<TH1F*> TPC_EnergyMC;  
  vector<TH1F*> TPC_EdepMC;
  vector<TH2F*> TPC_Veto_EnergyMC;  

  TNtuple *MCValue_ntuple;
  //  double c14_spectrum_endpoint;
  //  double c14_spectrum_startpoint;
  //  vector<double> c14_energy_spectrum, c14_energy_spectrum_1st, c14_energy_spectrum_2nd,c14_energy_spectrum_3rd;
  vector<double> fraction;
  //  TH1F **EnergyMC;
  //  vector<TNtuple*> VetoMC_ntuple;
    
  ClassDef(vetopulsesplit,0);
};

inline vector<int> vetopulsesplit::Colors()
{
  vector<int> color;
  color.push_back(TColor::GetColor("#FF2007")); //red         0
  color.push_back(TColor::GetColor("#5A1DE8")); //violet      1
  color.push_back(TColor::GetColor("#000000"));
  //  color.push_back(TColor::GetColor("#E8A60C")); //brown       2
  color.push_back(TColor::GetColor("#F73CFF")); //pink        5
  color.push_back(TColor::GetColor("#1CFFDF")); //low green   7
  //  color.push_back(TColor::GetColor("#FFFC19")); //yello       3
  color.push_back(TColor::GetColor("#1485CC")); //blue        4
  color.push_back(TColor::GetColor("#FF791F")); //orange      6
  color.push_back(TColor::GetColor("#AF2FCC")); //dark pink   8
  color.push_back(TColor::GetColor("#E8A60C"));
  color.push_back(TColor::GetColor("#B26618"));
  color.push_back(TColor::GetColor("#79FFFF"));
  color.push_back(TColor::GetColor("#11FF8F"));

  color.push_back(TColor::GetColor("#59FF49")); //green       12

  color.push_back(TColor::GetColor("#09B23E"));
  color.push_back(TColor::GetColor("#62B21D"));
  color.push_back(TColor::GetColor("#FB78FF"));
  color.push_back(TColor::GetColor("#AE3EFF"));
  color.push_back(TColor::GetColor("#B24F18"));
  color.push_back(TColor::GetColor("#FFA762"));
  color.push_back(TColor::GetColor("#2FDAFF"));
  color.push_back(TColor::GetColor("#1D537F"));
  color.push_back(TColor::GetColor("#4A067F"));
  
  return color;
}


#endif /*_vetopulsesplit_H*/
