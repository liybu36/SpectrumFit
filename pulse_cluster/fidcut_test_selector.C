#define fidcut_test_selector_cxx
// The class definition in fidcut_test_selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("fidcut_test_selector.C")
// Root > T->Process("fidcut_test_selector.C","some options")
// Root > T->Process("fidcut_test_selector.C+")
//

#include "fidcut_test_selector.h"
#include "TSystem.h"
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <cassert>
#include <sstream>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cstring>


void fidcut_test_selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   TString option = GetOption();

   xyl=0;


   fTree=0;
   fFile=0;
   fProofFile=0;


   
}

void fidcut_test_selector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
   xyl = (xylocator*) fInput->FindObject("xyl");
   TString option = GetOption();
  
   UInt_t opt = TProofOutputFile::kRegister | TProofOutputFile::kOverwrite | TProofOutputFile::kVerify;
   fProofFile = new TProofOutputFile(option,TProofOutputFile::kMerge, opt);

   fFile = fProofFile->OpenFile("RECREATE");
   if (fFile && fFile->IsZombie()) SafeDelete(fFile);

   fTree = new TTree("fidcut", "Tree with fiducial cut and other xy recon output");
   fTree->Branch("true_x",&true_x);
   fTree->Branch("true_y",&true_y);
   fTree->Branch("s2",&s2_copy);
   fTree->Branch("edge_chi2_ratio",&fidcut_ratio);
   fTree->Branch("best_x",&best_x);
   fTree->Branch("best_y",&best_y);
   fTree->Branch("recon_vec",&recon_vec);
   // File resident
   fTree->SetDirectory(fFile);
   fTree->AutoSave();

   evt_s2=0;
   true_x=0;
   true_y=0;

   recon_vec.clear();
   fidcut_ratio=0;
   best_x=-999;
   best_y=-999;



}

Bool_t fidcut_test_selector::Process(Long64_t entry)
{


   GetEntry(entry);

   double total_s2 = std::accumulate(evt_s2->begin(),evt_s2->end(),0.00, xyl->map_acc);

   recon_vec = xyl->leastsquares_vec(evt_s2,1);
   fidcut_ratio=xyl->edge_chi2_ratio(recon_vec);
   
   double best_chi2=1e99;
   for(std::vector<xychi2>::iterator vit = recon_vec.begin(); vit!= recon_vec.end(); ++vit){
     

      if(vit->chi2<best_chi2) {
         best_chi2 = vit->chi2;
         best_x=vit->x;
         best_y=vit->y;
      }
   }
   s2_copy=*evt_s2;
   fTree->Fill();
return kTRUE;
}

void fidcut_test_selector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   if (fFile) {
      if (!fTree) {
         Error("SlaveTerminate", "'tree' is undefined!");
         return;
      }
      Bool_t cleanup = kFALSE;
      TDirectory::TContext ctx(0);
      if (fTree->GetEntries() > 0) {
         fFile->cd();
         fTree->Write();
         fProofFile->Print();
         fOutput->Add(fProofFile);
      } else {
         cleanup = kTRUE;
      }
      fTree->SetDirectory(0);
      fFile->Close();
      // Cleanup, if needed
      if (cleanup) {
         TUrl uf(*(fFile->GetEndpointUrl()));
         SafeDelete(fFile);
         gSystem->Unlink(uf.GetFile());
         SafeDelete(fProofFile);
      }
   }
 }

void fidcut_test_selector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.      
}
