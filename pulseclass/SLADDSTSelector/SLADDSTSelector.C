#define SLADDSTSelector_cxx
// The class definition in SLADDSTSelector.h has been generated automatically
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
// Root > T->Process("SLADDSTSelector.C")
// Root > T->Process("SLADDSTSelector.C","some options")
// Root > T->Process("SLADDSTSelector.C+")
//

#include "SLADDSTSelector.h"
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
#include "TSystem.h"
using namespace std;

void SLADDSTSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void SLADDSTSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
  Info("SlaveBegin()","beginning .....");
  TString option = GetOption();

  string outdir="/darkside/users/hqian/pulse_splitter/VetoData/";
  string output = outdir+ option.Data();
  UInt_t opt = TProofOutputFile::kRegister | TProofOutputFile::kOverwrite | TProofOutputFile::kVerify;
  fProofFile = new TProofOutputFile(option,TProofOutputFile::kMerge, opt);
  //   fProofFile->SetOutputFileName(option.Data());
  fProofFile->SetOutputFileName(output.c_str());
  fFile = fProofFile->OpenFile("RECREATE");
  if (fFile && fFile->IsZombie()) SafeDelete(fFile);
  SetOutputTree(fFile);

  BookHistograms();
}

Bool_t SLADDSTSelector::Process(Long64_t entry)
{
  fChain->GetEntry(entry);
  FillHistograms();
  return kTRUE;
}

void SLADDSTSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
  if (fFile) {
    if (!lsvtree) {
      Error("SlaveTerminate", "'tree' is undefined!");
      return;
    }
    Bool_t cleanup = kFALSE;
    TDirectory *savedir = gDirectory;
    if (lsvtree->GetEntries() > 0) {
      fFile->cd();
      lsvtree->Write();
      fProofFile->Print();
      fOutput->Add(fProofFile);
    } else {
      cleanup = kTRUE;
    }
    lsvtree->SetDirectory(0);
    gDirectory = savedir;
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

void SLADDSTSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  TString option = GetOption();
  TList* list = GetOutputList();

  C14Spectrum = dynamic_cast<TH1F*>(list->FindObject(Form("C14Spectrum")));
  Co60Spectrum = dynamic_cast<TH1F*>(list->FindObject(Form("Co60Spectrum")));
  Fulllate = dynamic_cast<TH1F*>(list->FindObject(Form("Fulllate")));
  FullSpectrum = dynamic_cast<TH1F*>(list->FindObject(Form("FullSpectrum")));
  FullSpectrum_ntuple = dynamic_cast<TNtuple*>(list->FindObject(Form("FullSpectrum_ntuple")));
  Veto_TPC_hist = dynamic_cast<TH2F*>(list->FindObject(Form("Veto_TPC_hist")));
  id = dynamic_cast<TH1F*>(list->FindObject(Form("id")));
  int nEntries = id->Integral();

  TCanvas *c1 = new TCanvas("c1", "Full Spectrum with Fit", 1200, 600);
  c1->SetLogy();
  c1->cd();
  FullSpectrum->Draw();

  bool same = true;
  if(same){
    string outdir = "/darkside/users/hqian/pulse_splitter/VetoData/";
    string output =outdir+ option.Data();
    //    TFile outfile(output.c_str(), "RECREATE");
    TFile outfile(output.c_str(), "UPDATE");
    id->Write();
    C14Spectrum->Write();
    Co60Spectrum->Write();
    FullSpectrum->Write();
    Fulllate->Write();
    FullSpectrum_ntuple->Write();
    Veto_TPC_hist->Write();
    c1->Write();
    outfile.Write();
    outfile.Close();
  }
  else{
    if((fProofFile = dynamic_cast<TProofOutputFile*>(list->FindObject(option))))
      {
	//Write histograms into a root file
	//  TFile f(option, "RECREATE");
	TString outputFile( fProofFile->GetOutputFileName() );
	fFile = TFile::Open(outputFile);
	fFile->cd();
	id->Write();
	C14Spectrum->Write();
	Co60Spectrum->Write();
	FullSpectrum->Write();
	Fulllate->Write();
	FullSpectrum_ntuple->Write();
	Veto_TPC_hist->Write();
	c1->Write();
      }
    else {
      Error("Terminate", "TProofOutputFile not found");
      return;
    }
    fFile->Close();
  }

  Info("Terminate()","terminating ...%s",option.Data());

}

void SLADDSTSelector::SetOutputTree(TFile *fFile)
{
  lsvtree = new TTree("lsvtree","Reconstructed Events");
  lsvtree->Branch("charge", &charge,32000,1);
  //  lsvtree->Branch("multicut", &multicut,32000,1);
  lsvtree->Branch("width", &width,32000,1);
  lsvtree->Branch("height", &height,32000,1);
  //  lsvtree->Branch("pass", &pass,32000,1);

  lsvtree->SetDirectory(fFile);
  lsvtree->AutoSave();

}

void SLADDSTSelector::FillHistograms()
{
  id->Fill(0);
  //  for(int i=0; i<od_nclusters; i++)
  //TBAsymmetry correction
  Double_t par2[] = { -0.0397956, -0.27216, 0.794036,
		      1.70427, -3.98323, -8.50783, -2.66051 };
  Double_t TBAsym = (tpc_s1_top - tpc_s1_bottom) /
    (tpc_s1_top + tpc_s1_bottom);
  Double_t x = TBAsym; // normalized at 200 PE
  Double_t diff_total_s1 = par2[0] + (par2[1] + (par2[2] +
						 (par2[3] + (par2[4] + par2[5] * x) * x) * x) * x) *x;//(total_s1-total_s1_corr)/total_s1
  Double_t total_s1_TBAcorr = tpc_total_s1 * (1. - diff_total_s1);

  for(size_t i=0; i<od_cluster_start->size(); i++)
    {
      if(od_cluster_charge->at(i)>0 && od_cluster_pass_multcut->at(i)==1)
	{
	  double gps = od_cluster_dtprompt->at(i);
	  // if(TMath::Abs(gps)<0.05)
	  if(gps>-2 && gps<0.1)
	    {
	      FullSpectrum->Fill(od_cluster_charge->at(i));
	      C14Spectrum->Fill(od_cluster_charge->at(i));
	      FullSpectrum_ntuple->Fill(od_cluster_charge->at(i),od_cluster_height->at(i),od_cluster_multiplicity->at(i),gps);
	      charge=od_cluster_charge->at(i);
	      width=od_cluster_width->at(i);
	      height=od_cluster_height->at(i);
	      // multicut=od_cluster_multiplicity->at(i);
	      // pass=od_cluster_pass_multcut->at(i);
	      lsvtree->Fill();
	      Veto_TPC_hist->Fill(total_s1_TBAcorr,od_cluster_charge->at(i));
	      if(total_s1_TBAcorr>10200 && total_s1_TBAcorr<11600 && od_cluster_charge->at(i)>400 && od_cluster_charge->at(i)<800)
		Co60Spectrum->Fill(od_cluster_charge->at(i));
	    }
	  if(gps>40 && gps<140)
	    {
	      Fulllate->Fill(od_cluster_charge->at(i));
	    }
	}
    }
}

bool SLADDSTSelector::multicut(float height,float multiplicity, float charge){
  return height/multiplicity < (2.563e7 + TMath::Sqrt(1.574e14+1.390e12*TMath::Power(charge-14.40,2)));
}

void SLADDSTSelector::BookHistograms()
{
  Info("ReconSelector::BookHistograms()","creating histograms...");
  string label = GetOption();
  TList* list = GetOutputList();

  int Bins = 2000;

  Fulllate = new TH1F("Fulllate","lsv cluster fixed width charge full spectrum[40,140]us;charge [PE];Counts",Bins,0,2000);
  FullSpectrum = new TH1F("FullSpectrum","lsv cluster fixed width charge full spectrum;charge [PE];Counts",Bins,0,2000);
  FullSpectrum_ntuple = new TNtuple("FullSpectrum_ntuple","Veto All Data","charge:height:multi:time");
  C14Spectrum = new TH1F("C14Spectrum","lsv cluster fixed width charge;charge [PE];Counts",300,0,300);
  Co60Spectrum = new TH1F("Co60Spectrum","Co60 lsv cluster fixed width charge;charge [PE];Counts",1000,0,1000);
  Veto_TPC_hist = new TH2F("Veto_TPC_hist","lsv cluster charge vs tpc s1 TBAsym;total_s1_TBAcorr [PE];od_cluster_charge[PE]",1000,0,25.e+3,1000,0,2000);
  id = new TH1F("id","event ID",1,0,1);
  list->Add(id);
  list->Add(Fulllate);
  list->Add(FullSpectrum);
  list->Add(FullSpectrum_ntuple);
  list->Add(C14Spectrum);
  list->Add(Co60Spectrum);
  list->Add(Veto_TPC_hist);
}

