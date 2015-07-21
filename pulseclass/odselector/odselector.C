#define odselector_cxx
// The class definition in odselector.h has been generated automatically
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
// Root > T->Process("odselector.C")
// Root > T->Process("odselector.C","some options")
// Root > T->Process("odselector.C+")
//

#include "odselector.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <TH2.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TBox.h>
using namespace std;

//#define MoreBins

void odselector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void odselector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
   Info("SlaveBegin()","beginning .....");
   TString option = GetOption();
   BookHistograms();
   
}

Bool_t odselector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either odselector::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
  fChain->GetEntry(entry);
  FillHistograms();
  return kTRUE;
}

void odselector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void odselector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  string label = GetOption();
  TList* list = GetOutputList();
  C14Spectrum = dynamic_cast<TH1F*>(list->FindObject(Form("C14Spectrum")));
  FullSpectrum = dynamic_cast<TH1F*>(list->FindObject(Form("FullSpectrum")));
  FullSpectrum_ntuple = dynamic_cast<TNtuple*>(list->FindObject(Form("FullSpectrum_ntuple")));
  id = dynamic_cast<TH1F*>(list->FindObject(Form("id")));
  //  int nEntries = fChain->GetEntries();
  int nEntries = id->Integral();
  //  FullSpectrum->Sumw2();
  //  C14Spectrum->Sumw2();

  //  FullSpectrum->Scale(1.0/nEntries);
  //  C14Spectrum->Scale(1.0/nEntries);

  TCanvas *c1 = new TCanvas("c1", "Full Spectrum with Fit", 1200, 600);
  c1->SetLogy();
  c1->cd();
  FullSpectrum->Draw();

  //  string outdir = "/darkside/users/hqian/pulse_splitter/VetoData/":
  string output = label;
  /*
#ifdef MoreBins
  string output = outdir +"PulseSplitRealEnergy_MoreBins.root";
#else
  string output = outdir +"PulseSplitRealEnergy.root";
#endif
  TFile outfile(output.c_str(), "RECREATE");
  */
  TFile outfile(output.c_str(), "RECREATE");
  id->Write();
  C14Spectrum->Write();
  FullSpectrum->Write();
  FullSpectrum_ntuple->Write();
  c1->Write();
  outfile.Write();
  outfile.Close();

  Info("Terminate()","terminating ...%s",output.c_str());
}

void odselector::FillHistograms()
{
  id->Fill(0);
  for(int i=0; i<lsv_n_clusters; i++)
    {
      if(multicut(lsv_cluster_height->At(i),lsv_cluster_max_multiplicity->At(i),lsv_cluster_fixed_width_charge->At(i))
	 && lsv_cluster_start_ns->At(i)>3770 && lsv_cluster_start_ns->At(i)<3786 
	 && lsv_cluster_fixed_width_charge->At(i)>0 )
	{
	  FullSpectrum->Fill(lsv_cluster_fixed_width_charge->At(i));
	  FullSpectrum_ntuple->Fill(lsv_cluster_fixed_width_charge->At(i),lsv_cluster_height->At(i),lsv_cluster_max_multiplicity->At(i));
	  C14Spectrum->Fill(lsv_cluster_fixed_width_charge->At(i));
	}
    }
  
}

bool odselector::multicut(float height,float multiplicity, float charge){
  return height/multiplicity < (2.563e7 + TMath::Sqrt(1.574e14+1.390e12*TMath::Power(charge-14.40,2)));
}

void odselector::BookHistograms()
{
  Info("ReconSelector::BookHistograms()","creating histograms...");
  string label = GetOption();
  TList* list = GetOutputList();

#ifdef MoreBins
  int Bins = 2000;
#else
  int Bins = 1000;
#endif
  
  FullSpectrum = new TH1F("FullSpectrum","lsv cluster fixed width charge full spectrum;charge [PE];Counts",Bins,0,2000);
  FullSpectrum_ntuple = new TNtuple("FullSpectrum_ntuple","Veto All Data","charge:height:multi"); 
  C14Spectrum = new TH1F("C14Spectrum","lsv cluster fixed width charge;charge [PE];Counts",300,0,300);
  id = new TH1F("id","event ID",1,0,1);
  list->Add(id);
  list->Add(FullSpectrum);
  list->Add(FullSpectrum_ntuple);
  list->Add(C14Spectrum);
}
