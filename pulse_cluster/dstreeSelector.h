//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 25 13:37:41 2014 by ROOT version 5.34/21
// from TTree dstree/The G4DS Root Tree
// found on file: outtmb.root
//////////////////////////////////////////////////////////

#ifndef dstreeSelector_h
#define dstreeSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <vector>

#include "G4DSdata.hh"
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;

class dstreeSelector : public TSelector {
public :
  
  G4DSHitData currentHit;
  ReconEvent *recon_evts;
  std::vector<ReconEvent>* recon_evts;
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   // Declaration of leaf types
   Int_t           ev;
   Int_t           pdg;
   Float_t         ene0;
   Float_t         s1ene;
   Float_t         s2ene;
   Float_t         veto_visene;
   Float_t         mu_visene;
   Float_t         tpcene;
   Float_t         vetoene;
   Float_t         muene;
   Float_t         ene;
   Float_t         x;
   Float_t         y;
   Float_t         z;
   Float_t         r;
   Float_t         px;
   Float_t         py;
   Float_t         pz;
   Float_t         bx;
   Float_t         by;
   Float_t         bz;
   Int_t           npeS1;
   Int_t           npeS2;
   Int_t           npe;
   Int_t           munpe;
   Int_t           vnpe;
   Int_t           nph;
   Int_t           ndaughters;
   Int_t           ndeposits;
   Int_t           nusers;
   Int_t           dau_id[3];   //[ndaughters]
   Int_t           dau_pdg[3];   //[ndaughters]
   Int_t           dau_pid[3];   //[ndaughters]
   Int_t           dau_process[3];   //[ndaughters]
   Int_t           dau_mat[3];   //[ndaughters]
   Double_t        dau_time[3];   //[ndaughters]
   Float_t         dau_ene[3];   //[ndaughters]
   Float_t         dau_x[3];   //[ndaughters]
   Float_t         dau_y[3];   //[ndaughters]
   Float_t         dau_z[3];   //[ndaughters]
   Float_t         dau_r[3];   //[ndaughters]
   Float_t         dau_px[3];   //[ndaughters]
   Float_t         dau_py[3];   //[ndaughters]
   Float_t         dau_pz[3];   //[ndaughters]
   Int_t           dep_pdg[2831];   //[ndeposits]
   Int_t           dep_mat[2831];   //[ndeposits]
   Int_t           dep_track[2831];   //[ndeposits]
   Double_t        dep_time[2831];   //[ndeposits]
   Float_t         dep_ene[2831];   //[ndeposits]
   Float_t         dep_step[2831];   //[ndeposits]
   Float_t         dep_x[2831];   //[ndeposits]
   Float_t         dep_y[2831];   //[ndeposits]
   Float_t         dep_z[2831];   //[ndeposits]
   Float_t         dep_r[2831];   //[ndeposits]
   Int_t           userint1[1];   //[nusers]
   Int_t           userint2[1];   //[nusers]
   Float_t         userfloat1[1];   //[nusers]
   Float_t         userfloat2[1];   //[nusers]
   Double_t        userdouble0[1];   //[nusers]
   Double_t        pe_time[1];   //[npe]
   Int_t           pe_pmt[1];   //[npe]
   Double_t        vpe_time[1];   //[vnpe]
   Int_t           vpe_pmt[1];   //[vnpe]
   Double_t        mupe_time[1];   //[munpe]
   Int_t           mupe_pmt[1];   //[munpe]
   Int_t           ph_volume[1];   //[nph]
   Int_t           ph_pid[1];   //[nph]
   Float_t         ph_wl[1];   //[nph]
   Float_t         ph_x[1];   //[nph]
   Float_t         ph_y[1];   //[nph]
   Float_t         ph_z[1];   //[nph]
   Double_t        ph_time[1];   //[nph]

   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_ene0;   //!
   TBranch        *b_s1ene;   //!
   TBranch        *b_s2ene;   //!
   TBranch        *b_veto_visene;   //!
   TBranch        *b_mu_visene;   //!
   TBranch        *b_tpcene;   //!
   TBranch        *b_vetoene;   //!
   TBranch        *b_muene;   //!
   TBranch        *b_ene;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_radius;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_by;   //!
   TBranch        *b_bz;   //!
   TBranch        *b_npeS1;   //!
   TBranch        *b_npeS2;   //!
   TBranch        *b_npe;   //!
   TBranch        *b_munpe;   //!
   TBranch        *b_vnpe;   //!
   TBranch        *b_nph;   //!
   TBranch        *b_ndaughters;   //!
   TBranch        *b_ndeposits;   //!
   TBranch        *b_nusers;   //!
   TBranch        *b_dau_id;   //!
   TBranch        *b_dau_pdg;   //!
   TBranch        *b_dau_pid;   //!
   TBranch        *b_dau_process;   //!
   TBranch        *b_dau_mat;   //!
   TBranch        *b_dau_time;   //!
   TBranch        *b_dau_ene;   //!
   TBranch        *b_dau_x;   //!
   TBranch        *b_dau_y;   //!
   TBranch        *b_dau_z;   //!
   TBranch        *b_dau_r;   //!
   TBranch        *b_dau_px;   //!
   TBranch        *b_dau_py;   //!
   TBranch        *b_dau_pz;   //!
   TBranch        *b_dep_pdg;   //!
   TBranch        *b_dep_mat;   //!
   TBranch        *b_dep_track;   //!
   TBranch        *b_dep_time;   //!
   TBranch        *b_dep_ene;   //!
   TBranch        *b_dep_step;   //!
   TBranch        *b_dep_x;   //!
   TBranch        *b_dep_y;   //!
   TBranch        *b_dep_z;   //!
   TBranch        *b_dep_r;   //!
   TBranch        *b_userint1;   //!
   TBranch        *b_userint2;   //!
   TBranch        *b_userfloat1;   //!
   TBranch        *b_userfloat2;   //!
   TBranch        *b_userdouble0;   //!
   TBranch        *b_pe_time;   //!
   TBranch        *b_pe_pmt;   //!
   TBranch        *b_vpe_time;   //!
   TBranch        *b_vpe_pmt;   //!
   TBranch        *b_mupe_time;   //!
   TBranch        *b_mupe_pmt;   //!
   TBranch        *b_ph_volume;   //!
   TBranch        *b_ph_pid;   //!
   TBranch        *b_ph_wl;   //!
   TBranch        *b_ph_x;   //!
   TBranch        *b_ph_y;   //!
   TBranch        *b_ph_z;   //!
   TBranch        *b_ph_time;   //!

   int event_id, n_active;
   vector<double> *recon_et, *recon_ex, *recon_ey, *recon_ez, *recon_dt, *recon_dz, *recon_dr,
     *recon_edep, *recon_eqch, *recon_edep_electron, *recon_edep_nuclear, *recon_quenchingfactor, *recon_parent_enrg;
   vector<TString> *recon_vol, * recon_particle;
   vector<int> *recon_parent_pdg;
   vector<bool> *recon_contain_alpha;

 dstreeSelector(TTree * /*tree*/ =0) : fChain(0),
     recon_et(0), recon_ex(0), recon_ey(0), recon_ez(0), recon_dt(0), recon_dz(0), recon_dr(0),
     recon_edep(0), recon_eqch(0), recon_edep_electron(0), recon_edep_nuclear(0), recon_quenchingfactor(0), recon_parent_enrg(0),
     recon_vol(0), recon_particle(0), recon_parent_pdg(0), recon_contain_alpha(0), event_id(0), n_active(0), recon_evts(0);
     { }

   virtual ~dstreeSelector() { }
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
   TTree *ReconTree;

   void SetOutputTree();
   bool CompareRE(ReconEvent , ReconEvent);
   void InitializeRE(ReconEvent&,G4DSHitData&, int);
   bool CompatibleRE(ReconEvent&,ReconEvent);
   void MergeRE(ReconEvent&, ReconEvent);
   double get_point(vector<double>,vector<double>,double);
   double get_point(double*,double*,double,int);
   double o_scint_e_quenching(double);
   double get_quenching_factor_organic_scintillator(ReconEvent);
   double get_quenching_factor_liquid_argon(ReconEvent,bool);
   double get_quenching_factor(ReconEvent, bool);
   string convertPDG(int);
   bool checkActive(Float_t, Float_t);
   string convertVolume(int, Float_t, Float_t);
   string returnType(int);
   void FillHistograms(Long64_t);
   void returnHit(G4DSHitData&);
   //   int doRecon(const char*,string,string);
   ClassDef(dstreeSelector,0);
};

#endif

#ifdef dstreeSelector_cxx
void dstreeSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev, &b_ev);
   fChain->SetBranchAddress("pdg", &pdg, &b_pdg);
   fChain->SetBranchAddress("ene0", &ene0, &b_ene0);
   fChain->SetBranchAddress("s1ene", &s1ene, &b_s1ene);
   fChain->SetBranchAddress("s2ene", &s2ene, &b_s2ene);
   fChain->SetBranchAddress("veto_visene", &veto_visene, &b_veto_visene);
   fChain->SetBranchAddress("mu_visene", &mu_visene, &b_mu_visene);
   fChain->SetBranchAddress("tpcene", &tpcene, &b_tpcene);
   fChain->SetBranchAddress("vetoene", &vetoene, &b_vetoene);
   fChain->SetBranchAddress("muene", &muene, &b_muene);
   fChain->SetBranchAddress("ene", &ene, &b_ene);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("r", &r, &b_radius);
   fChain->SetBranchAddress("px", &px, &b_px);
   fChain->SetBranchAddress("py", &py, &b_py);
   fChain->SetBranchAddress("pz", &pz, &b_pz);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("by", &by, &b_by);
   fChain->SetBranchAddress("bz", &bz, &b_bz);
   fChain->SetBranchAddress("npeS1", &npeS1, &b_npeS1);
   fChain->SetBranchAddress("npeS2", &npeS2, &b_npeS2);
   fChain->SetBranchAddress("npe", &npe, &b_npe);
   fChain->SetBranchAddress("munpe", &munpe, &b_munpe);
   fChain->SetBranchAddress("vnpe", &vnpe, &b_vnpe);
   fChain->SetBranchAddress("nph", &nph, &b_nph);
   fChain->SetBranchAddress("ndaughters", &ndaughters, &b_ndaughters);
   fChain->SetBranchAddress("ndeposits", &ndeposits, &b_ndeposits);
   fChain->SetBranchAddress("nusers", &nusers, &b_nusers);
   fChain->SetBranchAddress("dau_id", dau_id, &b_dau_id);
   fChain->SetBranchAddress("dau_pdg", dau_pdg, &b_dau_pdg);
   fChain->SetBranchAddress("dau_pid", dau_pid, &b_dau_pid);
   fChain->SetBranchAddress("dau_process", dau_process, &b_dau_process);
   fChain->SetBranchAddress("dau_mat", dau_mat, &b_dau_mat);
   fChain->SetBranchAddress("dau_time", dau_time, &b_dau_time);
   fChain->SetBranchAddress("dau_ene", dau_ene, &b_dau_ene);
   fChain->SetBranchAddress("dau_x", dau_x, &b_dau_x);
   fChain->SetBranchAddress("dau_y", dau_y, &b_dau_y);
   fChain->SetBranchAddress("dau_z", dau_z, &b_dau_z);
   fChain->SetBranchAddress("dau_r", dau_r, &b_dau_r);
   fChain->SetBranchAddress("dau_px", dau_px, &b_dau_px);
   fChain->SetBranchAddress("dau_py", dau_py, &b_dau_py);
   fChain->SetBranchAddress("dau_pz", dau_pz, &b_dau_pz);
   fChain->SetBranchAddress("dep_pdg", dep_pdg, &b_dep_pdg);
   fChain->SetBranchAddress("dep_mat", dep_mat, &b_dep_mat);
   fChain->SetBranchAddress("dep_track", dep_track, &b_dep_track);
   fChain->SetBranchAddress("dep_time", dep_time, &b_dep_time);
   fChain->SetBranchAddress("dep_ene", dep_ene, &b_dep_ene);
   fChain->SetBranchAddress("dep_step", dep_step, &b_dep_step);
   fChain->SetBranchAddress("dep_x", dep_x, &b_dep_x);
   fChain->SetBranchAddress("dep_y", dep_y, &b_dep_y);
   fChain->SetBranchAddress("dep_z", dep_z, &b_dep_z);
   fChain->SetBranchAddress("dep_r", dep_r, &b_dep_r);
   fChain->SetBranchAddress("userint1", &userint1, &b_userint1);
   fChain->SetBranchAddress("userint2", &userint2, &b_userint2);
   fChain->SetBranchAddress("userfloat1", &userfloat1, &b_userfloat1);
   fChain->SetBranchAddress("userfloat2", &userfloat2, &b_userfloat2);
   fChain->SetBranchAddress("userdouble0", &userdouble0, &b_userdouble0);
   fChain->SetBranchAddress("pe_time", &pe_time, &b_pe_time);
   fChain->SetBranchAddress("pe_pmt", &pe_pmt, &b_pe_pmt);
   fChain->SetBranchAddress("vpe_time", &vpe_time, &b_vpe_time);
   fChain->SetBranchAddress("vpe_pmt", &vpe_pmt, &b_vpe_pmt);
   fChain->SetBranchAddress("mupe_time", &mupe_time, &b_mupe_time);
   fChain->SetBranchAddress("mupe_pmt", &mupe_pmt, &b_mupe_pmt);
   fChain->SetBranchAddress("ph_volume", &ph_volume, &b_ph_volume);
   fChain->SetBranchAddress("ph_pid", &ph_pid, &b_ph_pid);
   fChain->SetBranchAddress("ph_wl", &ph_wl, &b_ph_wl);
   fChain->SetBranchAddress("ph_x", &ph_x, &b_ph_x);
   fChain->SetBranchAddress("ph_y", &ph_y, &b_ph_y);
   fChain->SetBranchAddress("ph_z", &ph_z, &b_ph_z);
   fChain->SetBranchAddress("ph_time", &ph_time, &b_ph_time);
}

Bool_t dstreeSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef dstreeSelector_cxx
