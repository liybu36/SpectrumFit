//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 20 13:55:37 2015 by ROOT version 5.34/21
// from TTree odtree/DarkSide-50 Outer Detector Event Tree
// found on file: ODRun000177.root
//////////////////////////////////////////////////////////

#ifndef odselector_h
#define odselector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TNtuple.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <TArrayF.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;
class odselector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event_number;
   UShort_t        trg_id;
   UShort_t        pps_counter;
   UInt_t          gps_fine_time_counter;
   UInt_t          gps_one_second_counter;
   UInt_t          bad_time_alignment;
   Double_t        absolute_t[64];
   Double_t        relative_t[64];
   Float_t         peak_time[256];
   Int_t           raw_npmts;
   Int_t           raw_npulses;
   Int_t           ch_raw_npulses[256];
   Float_t         lsv_total_spe_charge;
   Float_t         wt_total_spe_charge;
   Float_t         chan_spe_charge[256];
   vector<vector<double> > *all_spe_charges;
   Float_t         lsv_total_cluster_charge;
   Float_t         wt_total_cluster_charge;
   Float_t         lsv_max_cluster_charge;
   Float_t         wt_max_cluster_charge;
   Float_t         lsv_total_cluster_fixed_width_charge;
   Float_t         wt_total_cluster_fixed_width_charge;
   Int_t           lsv_n_clusters;
   Int_t           wt_n_clusters;
   TArrayF         *lsv_cluster_charge;
   TArrayF         *lsv_cluster_fixed_width_charge;
   TArrayF         *lsv_cluster_sat_corr_charge;
   TArrayF         *lsv_cluster_start_ns;
   TArrayF         *lsv_cluster_end_ns;
   TArrayF         *lsv_cluster_width_ns;
   TArrayF         *lsv_cluster_peak_ns;
   TArrayF         *lsv_cluster_peak_multiplicity;
   TArrayF         *lsv_cluster_max_multiplicity;
   TArrayF         *lsv_cluster_height;
   TArrayF         *lsv_cluster_frac_above;
   TArrayF         *lsv_cluster_saturated_chan;
   TArrayF         *lsv_cluster_dispersion;
   TArrayF         *lsv_cluster_theta;
   TArrayF         *lsv_cluster_phi;
   TArrayF         *lsv_cluster_x;
   TArrayF         *lsv_cluster_y;
   TArrayF         *lsv_cluster_z;
   TArrayF         *lsv_cluster_r;
   TArrayF         *lsv_cluster_early_integral;
   TArrayF         *wt_cluster_charge;
   TArrayF         *wt_cluster_fixed_width_charge;
   TArrayF         *wt_cluster_start_ns;
   TArrayF         *wt_cluster_end_ns;
   TArrayF         *wt_cluster_width_ns;
   TArrayF         *wt_cluster_peak_ns;
   TArrayF         *wt_cluster_peak_multiplicity;
   TArrayF         *wt_cluster_max_multiplicity;
   TArrayF         *wt_cluster_height;
   TArrayF         *wt_cluster_saturated_chan;
   TArrayF         *wt_cluster_dispersion;
   Float_t         lsv_max_window_charge;
   Float_t         lsv_max_window_height;
   Float_t         lsv_max_window_time;
   Float_t         wt_max_window_charge;
   Float_t         wt_max_window_height;
   Float_t         wt_max_window_time;
   Float_t         lsv_delayed_window_charge;
   Float_t         lsv_delayed_window_height;
   Float_t         lsv_delayed_window_time;
   Float_t         wt_delayed_window_charge;
   Float_t         wt_delayed_window_height;
   Float_t         wt_delayed_window_time;
   Float_t         lsv_late_window_charge;
   Float_t         lsv_late_window_height;
   Float_t         lsv_late_window_time;
   Float_t         wt_late_window_charge;
   Float_t         wt_late_window_height;
   Float_t         wt_late_window_time;
   Float_t         lsv_max_multiplicity;
   Float_t         lsv_max_multiplicity_time;
   Float_t         lsv_delayed_multiplicity;
   Float_t         lsv_delayed_multiplicity_time;
   Float_t         lsv_late_multiplicity;
   Float_t         lsv_late_multiplicity_time;
   Float_t         wt_max_multiplicity;
   Float_t         wt_max_multiplicity_time;
   Float_t         wt_delayed_multiplicity;
   Float_t         wt_delayed_multiplicity_time;
   Float_t         wt_late_multiplicity;
   Float_t         wt_late_multiplicity_time;
   Float_t         roi_charge[2];
   Float_t         roi_multiplicity[2];
   Float_t         roi_height[2];

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event_number;   //!
   TBranch        *b_trg_id;   //!
   TBranch        *b_pps_counter;   //!
   TBranch        *b_gps_fine_time_counter;   //!
   TBranch        *b_gps_one_second_counter;   //!
   TBranch        *b_bad_time_alignment;   //!
   TBranch        *b_absolute_t;   //!
   TBranch        *b_relative_t;   //!
   TBranch        *b_peak_time;   //!
   TBranch        *b_raw_npmts;   //!
   TBranch        *b_raw_npulses;   //!
   TBranch        *b_ch_raw_npulses;   //!
   TBranch        *b_lsv_total_spe_charge;   //!
   TBranch        *b_wt_total_spe_charge;   //!
   TBranch        *b_chan_spe_charge;   //!
   TBranch        *b_all_spe_charges;   //!
   TBranch        *b_lsv_total_cluster_charge;   //!
   TBranch        *b_wt_total_cluster_charge;   //!
   TBranch        *b_lsv_max_cluster_charge;   //!
   TBranch        *b_wt_max_cluster_charge;   //!
   TBranch        *b_lsv_total_cluster_fixed_width_charge;   //!
   TBranch        *b_wt_total_cluster_fixed_width_charge;   //!
   TBranch        *b_lsv_n_clusters;   //!
   TBranch        *b_wt_n_clusters;   //!
   TBranch        *b_lsv_cluster_charge;   //!
   TBranch        *b_lsv_cluster_fixed_width_charge;   //!
   TBranch        *b_lsv_cluster_sat_corr_charge;   //!
   TBranch        *b_lsv_cluster_start_ns;   //!
   TBranch        *b_lsv_cluster_end_ns;   //!
   TBranch        *b_lsv_cluster_width_ns;   //!
   TBranch        *b_lsv_cluster_peak_ns;   //!
   TBranch        *b_lsv_cluster_peak_multiplicity;   //!
   TBranch        *b_lsv_cluster_max_multiplicity;   //!
   TBranch        *b_lsv_cluster_height;   //!
   TBranch        *b_lsv_cluster_frac_above;   //!
   TBranch        *b_lsv_cluster_saturated_chan;   //!
   TBranch        *b_lsv_cluster_dispersion;   //!
   TBranch        *b_lsv_cluster_theta;   //!
   TBranch        *b_lsv_cluster_phi;   //!
   TBranch        *b_lsv_cluster_x;   //!
   TBranch        *b_lsv_cluster_y;   //!
   TBranch        *b_lsv_cluster_z;   //!
   TBranch        *b_lsv_cluster_r;   //!
   TBranch        *b_lsv_cluster_early_integral;   //!
   TBranch        *b_wt_cluster_charge;   //!
   TBranch        *b_wt_cluster_fixed_width_charge;   //!
   TBranch        *b_wt_cluster_start_ns;   //!
   TBranch        *b_wt_cluster_end_ns;   //!
   TBranch        *b_wt_cluster_width_ns;   //!
   TBranch        *b_wt_cluster_peak_ns;   //!
   TBranch        *b_wt_cluster_peak_multiplicity;   //!
   TBranch        *b_wt_cluster_max_multiplicity;   //!
   TBranch        *b_wt_cluster_height;   //!
   TBranch        *b_wt_cluster_saturated_chan;   //!
   TBranch        *b_wt_cluster_dispersion;   //!
   TBranch        *b_lsv_max_window_charge;   //!
   TBranch        *b_lsv_max_window_height;   //!
   TBranch        *b_lsv_max_window_time;   //!
   TBranch        *b_wt_max_window_charge;   //!
   TBranch        *b_wt_max_window_height;   //!
   TBranch        *b_wt_max_window_time;   //!
   TBranch        *b_lsv_delayed_window_charge;   //!
   TBranch        *b_lsv_delayed_window_height;   //!
   TBranch        *b_lsv_delayed_window_time;   //!
   TBranch        *b_wt_delayed_window_charge;   //!
   TBranch        *b_wt_delayed_window_height;   //!
   TBranch        *b_wt_delayed_window_time;   //!
   TBranch        *b_lsv_late_window_charge;   //!
   TBranch        *b_lsv_late_window_height;   //!
   TBranch        *b_lsv_late_window_time;   //!
   TBranch        *b_wt_late_window_charge;   //!
   TBranch        *b_wt_late_window_height;   //!
   TBranch        *b_wt_late_window_time;   //!
   TBranch        *b_lsv_max_multiplicity;   //!
   TBranch        *b_lsv_max_multiplicity_time;   //!
   TBranch        *b_lsv_delayed_multiplicity;   //!
   TBranch        *b_lsv_delayed_multiplicity_time;   //!
   TBranch        *b_lsv_late_multiplicity;   //!
   TBranch        *b_lsv_late_multiplicity_time;   //!
   TBranch        *b_wt_max_multiplicity;   //!
   TBranch        *b_wt_max_multiplicity_time;   //!
   TBranch        *b_wt_delayed_multiplicity;   //!
   TBranch        *b_wt_delayed_multiplicity_time;   //!
   TBranch        *b_wt_late_multiplicity;   //!
   TBranch        *b_wt_late_multiplicity_time;   //!
   TBranch        *b_roi_charge;   //!
   TBranch        *b_roi_multiplicity;   //!
   TBranch        *b_roi_height;   //!

 odselector(TTree * /*tree*/ =0) : fChain(0),id(0),
     FullSpectrum(0),FullSpectrum_ntuple(0)
     { }
   virtual ~odselector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

 private:
   void BookHistograms();
   void FillHistograms();
   bool multicut(float,float,float);
   TH1F *id;
   TH1F *FullSpectrum;
   TNtuple *FullSpectrum_ntuple;

   ClassDef(odselector,0);
};

#endif

#ifdef odselector_cxx
void odselector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   all_spe_charges = 0;
   lsv_cluster_charge = 0;
   lsv_cluster_fixed_width_charge = 0;
   lsv_cluster_sat_corr_charge = 0;
   lsv_cluster_start_ns = 0;
   lsv_cluster_end_ns = 0;
   lsv_cluster_width_ns = 0;
   lsv_cluster_peak_ns = 0;
   lsv_cluster_peak_multiplicity = 0;
   lsv_cluster_max_multiplicity = 0;
   lsv_cluster_height = 0;
   lsv_cluster_frac_above = 0;
   lsv_cluster_saturated_chan = 0;
   lsv_cluster_dispersion = 0;
   lsv_cluster_theta = 0;
   lsv_cluster_phi = 0;
   lsv_cluster_x = 0;
   lsv_cluster_y = 0;
   lsv_cluster_z = 0;
   lsv_cluster_r = 0;
   lsv_cluster_early_integral = 0;
   wt_cluster_charge = 0;
   wt_cluster_fixed_width_charge = 0;
   wt_cluster_start_ns = 0;
   wt_cluster_end_ns = 0;
   wt_cluster_width_ns = 0;
   wt_cluster_peak_ns = 0;
   wt_cluster_peak_multiplicity = 0;
   wt_cluster_max_multiplicity = 0;
   wt_cluster_height = 0;
   wt_cluster_saturated_chan = 0;
   wt_cluster_dispersion = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event_number", &event_number, &b_event_number);
   fChain->SetBranchAddress("trg_id", &trg_id, &b_trg_id);
   fChain->SetBranchAddress("pps_counter", &pps_counter, &b_pps_counter);
   fChain->SetBranchAddress("gps_fine_time_counter", &gps_fine_time_counter, &b_gps_fine_time_counter);
   fChain->SetBranchAddress("gps_one_second_counter", &gps_one_second_counter, &b_gps_one_second_counter);
   fChain->SetBranchAddress("bad_time_alignment", &bad_time_alignment, &b_bad_time_alignment);
   fChain->SetBranchAddress("absolute_t", absolute_t, &b_absolute_t);
   fChain->SetBranchAddress("relative_t", relative_t, &b_relative_t);
   fChain->SetBranchAddress("peak_time", peak_time, &b_peak_time);
   fChain->SetBranchAddress("raw_npmts", &raw_npmts, &b_raw_npmts);
   fChain->SetBranchAddress("raw_npulses", &raw_npulses, &b_raw_npulses);
   fChain->SetBranchAddress("ch_raw_npulses", ch_raw_npulses, &b_ch_raw_npulses);
   fChain->SetBranchAddress("lsv_total_spe_charge", &lsv_total_spe_charge, &b_lsv_total_spe_charge);
   fChain->SetBranchAddress("wt_total_spe_charge", &wt_total_spe_charge, &b_wt_total_spe_charge);
   fChain->SetBranchAddress("chan_spe_charge", chan_spe_charge, &b_chan_spe_charge);
   fChain->SetBranchAddress("all_spe_charges", &all_spe_charges, &b_all_spe_charges);
   fChain->SetBranchAddress("lsv_total_cluster_charge", &lsv_total_cluster_charge, &b_lsv_total_cluster_charge);
   fChain->SetBranchAddress("wt_total_cluster_charge", &wt_total_cluster_charge, &b_wt_total_cluster_charge);
   fChain->SetBranchAddress("lsv_max_cluster_charge", &lsv_max_cluster_charge, &b_lsv_max_cluster_charge);
   fChain->SetBranchAddress("wt_max_cluster_charge", &wt_max_cluster_charge, &b_wt_max_cluster_charge);
   fChain->SetBranchAddress("lsv_total_cluster_fixed_width_charge", &lsv_total_cluster_fixed_width_charge, &b_lsv_total_cluster_fixed_width_charge);
   fChain->SetBranchAddress("wt_total_cluster_fixed_width_charge", &wt_total_cluster_fixed_width_charge, &b_wt_total_cluster_fixed_width_charge);
   fChain->SetBranchAddress("lsv_n_clusters", &lsv_n_clusters, &b_lsv_n_clusters);
   fChain->SetBranchAddress("wt_n_clusters", &wt_n_clusters, &b_wt_n_clusters);
   fChain->SetBranchAddress("lsv_cluster_charge", &lsv_cluster_charge, &b_lsv_cluster_charge);
   fChain->SetBranchAddress("lsv_cluster_fixed_width_charge", &lsv_cluster_fixed_width_charge, &b_lsv_cluster_fixed_width_charge);
   fChain->SetBranchAddress("lsv_cluster_sat_corr_charge", &lsv_cluster_sat_corr_charge, &b_lsv_cluster_sat_corr_charge);
   fChain->SetBranchAddress("lsv_cluster_start_ns", &lsv_cluster_start_ns, &b_lsv_cluster_start_ns);
   fChain->SetBranchAddress("lsv_cluster_end_ns", &lsv_cluster_end_ns, &b_lsv_cluster_end_ns);
   fChain->SetBranchAddress("lsv_cluster_width_ns", &lsv_cluster_width_ns, &b_lsv_cluster_width_ns);
   fChain->SetBranchAddress("lsv_cluster_peak_ns", &lsv_cluster_peak_ns, &b_lsv_cluster_peak_ns);
   fChain->SetBranchAddress("lsv_cluster_peak_multiplicity", &lsv_cluster_peak_multiplicity, &b_lsv_cluster_peak_multiplicity);
   fChain->SetBranchAddress("lsv_cluster_max_multiplicity", &lsv_cluster_max_multiplicity, &b_lsv_cluster_max_multiplicity);
   fChain->SetBranchAddress("lsv_cluster_height", &lsv_cluster_height, &b_lsv_cluster_height);
   fChain->SetBranchAddress("lsv_cluster_frac_above", &lsv_cluster_frac_above, &b_lsv_cluster_frac_above);
   fChain->SetBranchAddress("lsv_cluster_saturated_chan", &lsv_cluster_saturated_chan, &b_lsv_cluster_saturated_chan);
   fChain->SetBranchAddress("lsv_cluster_dispersion", &lsv_cluster_dispersion, &b_lsv_cluster_dispersion);
   fChain->SetBranchAddress("lsv_cluster_theta", &lsv_cluster_theta, &b_lsv_cluster_theta);
   fChain->SetBranchAddress("lsv_cluster_phi", &lsv_cluster_phi, &b_lsv_cluster_phi);
   fChain->SetBranchAddress("lsv_cluster_x", &lsv_cluster_x, &b_lsv_cluster_x);
   fChain->SetBranchAddress("lsv_cluster_y", &lsv_cluster_y, &b_lsv_cluster_y);
   fChain->SetBranchAddress("lsv_cluster_z", &lsv_cluster_z, &b_lsv_cluster_z);
   fChain->SetBranchAddress("lsv_cluster_r", &lsv_cluster_r, &b_lsv_cluster_r);
   fChain->SetBranchAddress("lsv_cluster_early_integral", &lsv_cluster_early_integral, &b_lsv_cluster_early_integral);
   fChain->SetBranchAddress("wt_cluster_charge", &wt_cluster_charge, &b_wt_cluster_charge);
   fChain->SetBranchAddress("wt_cluster_fixed_width_charge", &wt_cluster_fixed_width_charge, &b_wt_cluster_fixed_width_charge);
   fChain->SetBranchAddress("wt_cluster_start_ns", &wt_cluster_start_ns, &b_wt_cluster_start_ns);
   fChain->SetBranchAddress("wt_cluster_end_ns", &wt_cluster_end_ns, &b_wt_cluster_end_ns);
   fChain->SetBranchAddress("wt_cluster_width_ns", &wt_cluster_width_ns, &b_wt_cluster_width_ns);
   fChain->SetBranchAddress("wt_cluster_peak_ns", &wt_cluster_peak_ns, &b_wt_cluster_peak_ns);
   fChain->SetBranchAddress("wt_cluster_peak_multiplicity", &wt_cluster_peak_multiplicity, &b_wt_cluster_peak_multiplicity);
   fChain->SetBranchAddress("wt_cluster_max_multiplicity", &wt_cluster_max_multiplicity, &b_wt_cluster_max_multiplicity);
   fChain->SetBranchAddress("wt_cluster_height", &wt_cluster_height, &b_wt_cluster_height);
   fChain->SetBranchAddress("wt_cluster_saturated_chan", &wt_cluster_saturated_chan, &b_wt_cluster_saturated_chan);
   fChain->SetBranchAddress("wt_cluster_dispersion", &wt_cluster_dispersion, &b_wt_cluster_dispersion);
   fChain->SetBranchAddress("lsv_max_window_charge", &lsv_max_window_charge, &b_lsv_max_window_charge);
   fChain->SetBranchAddress("lsv_max_window_height", &lsv_max_window_height, &b_lsv_max_window_height);
   fChain->SetBranchAddress("lsv_max_window_time", &lsv_max_window_time, &b_lsv_max_window_time);
   fChain->SetBranchAddress("wt_max_window_charge", &wt_max_window_charge, &b_wt_max_window_charge);
   fChain->SetBranchAddress("wt_max_window_height", &wt_max_window_height, &b_wt_max_window_height);
   fChain->SetBranchAddress("wt_max_window_time", &wt_max_window_time, &b_wt_max_window_time);
   fChain->SetBranchAddress("lsv_delayed_window_charge", &lsv_delayed_window_charge, &b_lsv_delayed_window_charge);
   fChain->SetBranchAddress("lsv_delayed_window_height", &lsv_delayed_window_height, &b_lsv_delayed_window_height);
   fChain->SetBranchAddress("lsv_delayed_window_time", &lsv_delayed_window_time, &b_lsv_delayed_window_time);
   fChain->SetBranchAddress("wt_delayed_window_charge", &wt_delayed_window_charge, &b_wt_delayed_window_charge);
   fChain->SetBranchAddress("wt_delayed_window_height", &wt_delayed_window_height, &b_wt_delayed_window_height);
   fChain->SetBranchAddress("wt_delayed_window_time", &wt_delayed_window_time, &b_wt_delayed_window_time);
   fChain->SetBranchAddress("lsv_late_window_charge", &lsv_late_window_charge, &b_lsv_late_window_charge);
   fChain->SetBranchAddress("lsv_late_window_height", &lsv_late_window_height, &b_lsv_late_window_height);
   fChain->SetBranchAddress("lsv_late_window_time", &lsv_late_window_time, &b_lsv_late_window_time);
   fChain->SetBranchAddress("wt_late_window_charge", &wt_late_window_charge, &b_wt_late_window_charge);
   fChain->SetBranchAddress("wt_late_window_height", &wt_late_window_height, &b_wt_late_window_height);
   fChain->SetBranchAddress("wt_late_window_time", &wt_late_window_time, &b_wt_late_window_time);
   fChain->SetBranchAddress("lsv_max_multiplicity", &lsv_max_multiplicity, &b_lsv_max_multiplicity);
   fChain->SetBranchAddress("lsv_max_multiplicity_time", &lsv_max_multiplicity_time, &b_lsv_max_multiplicity_time);
   fChain->SetBranchAddress("lsv_delayed_multiplicity", &lsv_delayed_multiplicity, &b_lsv_delayed_multiplicity);
   fChain->SetBranchAddress("lsv_delayed_multiplicity_time", &lsv_delayed_multiplicity_time, &b_lsv_delayed_multiplicity_time);
   fChain->SetBranchAddress("lsv_late_multiplicity", &lsv_late_multiplicity, &b_lsv_late_multiplicity);
   fChain->SetBranchAddress("lsv_late_multiplicity_time", &lsv_late_multiplicity_time, &b_lsv_late_multiplicity_time);
   fChain->SetBranchAddress("wt_max_multiplicity", &wt_max_multiplicity, &b_wt_max_multiplicity);
   fChain->SetBranchAddress("wt_max_multiplicity_time", &wt_max_multiplicity_time, &b_wt_max_multiplicity_time);
   fChain->SetBranchAddress("wt_delayed_multiplicity", &wt_delayed_multiplicity, &b_wt_delayed_multiplicity);
   fChain->SetBranchAddress("wt_delayed_multiplicity_time", &wt_delayed_multiplicity_time, &b_wt_delayed_multiplicity_time);
   fChain->SetBranchAddress("wt_late_multiplicity", &wt_late_multiplicity, &b_wt_late_multiplicity);
   fChain->SetBranchAddress("wt_late_multiplicity_time", &wt_late_multiplicity_time, &b_wt_late_multiplicity_time);
   fChain->SetBranchAddress("roi_charge", roi_charge, &b_roi_charge);
   fChain->SetBranchAddress("roi_multiplicity", roi_multiplicity, &b_roi_multiplicity);
   fChain->SetBranchAddress("roi_height", roi_height, &b_roi_height);
}

Bool_t odselector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef odselector_cxx
