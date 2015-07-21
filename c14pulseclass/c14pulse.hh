#ifndef _c14pulse_H
#define _c14pulse_H

#include "TH1F.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TRint.h"
#include "TColor.h"

using namespace std;
class c14pulse:public TObject 
{
public :
  c14pulse():
    C14_EnergyMC_1st(0),C14_EnergyMC_2nd(0), C14_EnergyMC_3rd(0),C14_Quenching_Model(0),
    quenching_plot(0)
  { c14_spectrum_endpoint = 156.27; 
    c14_spectrum_startpoint = 0;     
  }
  virtual ~c14pulse() { 
    MCinputdir.clear();
    MCDataFile.clear();
    RealDataFile.clear();
    Realinputdir.clear();
    outputdir.clear();
    c14_energy_spectrum.clear();
    c14_energy_spectrum_1st.clear();
    c14_energy_spectrum_2nd.clear();
    c14_energy_spectrum_3rd.clear();
    C14_EnergyMC_1st=NULL;
    C14_EnergyMC_2nd=NULL;
    C14_EnergyMC_3rd=NULL;
    C14_Quenching_Model=NULL;
    quenching_plot=NULL;
    BGFraction.clear();
  }
  void Recon_DataAnalysis();
  void ODTree_DataAnalysis(int,int);
  void Calculate_MCData(TTree*,TH1F*);
  void Readdatafile(TChain*,int);
  void Readdatafile(TChain*, int, int);
  void Calculate_Convolution(TH1F*, TH1F*, TH1F*);
  void Calculate_Norm(TH1F*, vector<double> &);
  double scint_e_quenching(double, double*);
  double scint_e_quenching(double);
  void Scint_Quenching_Model(TTree*,TH1F*,TH1F*);
  double response_function(double, double, double*,int);
  // int Main_Fit(int, int,string);
  double C14Fit(double*,double*);
  double C14Fit_1st(double*,double*);
  double C14Fit_2nd(double*,double*);

  /*  double Co60Fit(double*,double*);
  double Co57Fit(double*,double*);
  double K40Fit(double*,double*);
  double Tl208Fit(double*,double*);
  double Th232Fit(double*,double*);
  double U238Fit(double*,double*);  
  */
  double TotalFit(double*,double*);  
  TObjArray GetHistArray() { return Histlist; }
  string GetMCinputdir() { return MCinputdir; }
  void SetMCinputdir(string val) { MCinputdir=val; }
  string GetMCDataFile() { return MCDataFile; }
  void SetMCDataFile(string val) { MCDataFile=val; }

  string GetRealinputdir() { return Realinputdir; }
  void SetRealinputdir(string val) { Realinputdir=val; }
  string GetRealDataFile() { return RealDataFile; }
  void SetRealDataFile(string val) { RealDataFile=val; }

  string Getoutputdir() { return outputdir; }
  void Setoutputdir(string val) { outputdir=val; }
  vector<int> Colors();
  double Isotope_Sum(int,double*);
  void Isotope_Fraction(int,double*,double*);
  bool multicut(float,float,float);
  vector<double> BGFraction;

private:
  string MCinputdir;
  string MCDataFile;
  string Realinputdir;
  string RealDataFile;
  string outputdir;
  double c14_spectrum_endpoint;
  double c14_spectrum_startpoint;
  vector<double> c14_energy_spectrum, c14_energy_spectrum_1st, c14_energy_spectrum_2nd,c14_energy_spectrum_3rd;
  TH1F **EnergyMC;
  TH1F *C14_EnergyMC_1st, *C14_EnergyMC_2nd, *C14_EnergyMC_3rd;
  TH1F *C14_Quenching_Model, *quenching_plot;
  //  TH1F *C14_EnergyMC,*Co60_EnergyMC,*Co57_EnergyMC,*Tl208_EnergyMC,*K40_EnergyMC,*Th232_EnergyMC,*U238_EnergyMC;
  TObjArray Histlist;
  //  TH1F *Full_Spectrum;
  
  //  TF1 *Fit_C14, *Fit_C14_1st, *Fit_C14_2nd, *Fit_C14_3rd;
  //  TF1 *Fit_Co60,*Fit_Co57,*Fit_Tl208,*Fit_K40,*Fit_Th232,*Fit_U238;
  //  TF1 *Fit_Total;  
  
  ClassDef(c14pulse,0);
};

inline vector<int> c14pulse::Colors()
{
  vector<int> color;
  color.push_back(TColor::GetColor("#FF0000")); //red
  color.push_back(TColor::GetColor("#5A1DE8")); //violet
  color.push_back(TColor::GetColor("#1CFFDF")); //low green
  color.push_back(TColor::GetColor("#E8A60C")); //brown
  color.push_back(TColor::GetColor("#59FF49")); //green
  color.push_back(TColor::GetColor("#FFFC19")); //yello
  color.push_back(TColor::GetColor("#1485CC")); //blue
  color.push_back(TColor::GetColor("#F73CFF")); //pink
  color.push_back(TColor::GetColor("#FF791F")); //orange
  color.push_back(TColor::GetColor("#AF2FCC")); //dark pink

  color.push_back(TColor::GetColor("#E8A60C"));
  color.push_back(TColor::GetColor("#FF0000"));
  color.push_back(TColor::GetColor("#5A1DE8"));
  color.push_back(TColor::GetColor("#1CFFDF"));
  color.push_back(TColor::GetColor("#B26618"));
  color.push_back(TColor::GetColor("#FF379F"));
  color.push_back(TColor::GetColor("#25DCFF"));
  color.push_back(TColor::GetColor("#92B20B"));
  color.push_back(TColor::GetColor("#79FFFF"));
  color.push_back(TColor::GetColor("#11FF8F"));
  color.push_back(TColor::GetColor("#09B23E"));
  color.push_back(TColor::GetColor("#62B21D"));
  color.push_back(TColor::GetColor("#FB78FF"));
  color.push_back(TColor::GetColor("#AE3EFF"));
  color.push_back(TColor::GetColor("#B24F18"));
  color.push_back(TColor::GetColor("#FFA762"));
  color.push_back(TColor::GetColor("#2FDAFF"));
  color.push_back(TColor::GetColor("#1D537F"));
  color.push_back(TColor::GetColor("#4A067F"));
  color.push_back(TColor::GetColor("#000000"));
  
  return color;
}


#endif /*_c14pulse_H*/
