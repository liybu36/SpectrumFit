//It's the clustering algorithm using the TSelector.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TRint.h"
#include "TProof.h"
#include "TDSet.h"

using namespace std;
#include "reconSelector.h"

//const string Time = "_Oct13AM";
TRint* theApp;

//string dirname="/darkside/users/hqian/TMBfraction/";
//string dirname="/darkside/users/hqian/VetoEfficiency/";
string dirname="/darkside/users/hqian/60KeVGammaStudy/";

void Readdatafile(TChain *t, int start, int end)
{ 
  //  string dirname="/darkside/users/hqian/TMBfraction/";
  string filename;
  string middle="outtmb";
  string last=".root";
  stringstream oss;
  for(int i=start; i<=end; i++)
    {
      if(i==0)
	filename=dirname+middle+last;
      else
	{
	  oss<<i;
	  filename=dirname+middle+"_v"+oss.str()+last;
	}
      t->Add(filename.c_str());
      cout<<"Processing Data file: "<<filename<<endl;
      oss.str("");      
    }
}

#define test_data
string Readdatafile(int i)
{ 
  //  string dirname="/cache/shared/darkside/hqian/TMBfraction/";
  string filename;
#ifdef test_data
  //  string middle="outneutron";
  string middle="out60gamma";
#else
  string middle="outtmb";
#endif
  string last=".root";
  stringstream oss;
  if(i==0)
    filename=dirname+middle+last;
  else
    {
      oss<<i;
      filename=dirname+middle+"_v"+oss.str()+last;
    }
  cout<<"Processing Data file: "<<filename<<endl;
  //  TFile *file = new TFile(filename.c_str());
  // t = (TTree*)file->Get("dstree");
  oss.str("");
  return filename;
}

void Process(TChain* chain, string label){
  TString option = label;
  //  cout<<"option: "<<option<<" \t label"<<label<<endl;
#define use_TProof
#ifdef use_TProof
  chain->SetProof();
  TProof* pr = TProof::Open("workers=2");
  //  pr->SetLogLevel(2, TProofDebug::kPacketizer);
  pr->SetParameter("PROOF_Packetizer","TPacketizer");                                                                             
  pr->SetParameter("PROOF_MaxSlavesPerNode",8);       
  chain->Process("reconSelector.C+",option.Data(),-1,0);  
#else
  reconSelector *selector = new reconSelector();
  selector->SetTRint(theApp);
  chain->Process(selector,option.Data());
  gSystem->Exit(0);
#endif
  //Print more information                                                                                                           
  //  pr->SetLogLevel(2,TProofDebug::kPacketizer);                                                                                   
  // pr->SetParameter("PROOF_Packetizer","TPacketizer");                                                                             
  // pr->SetParameter("PROOF_MaxSlavesPerNode",8);       
  //  pr->Mgr("palpatine.Princeton.EDU")->GetSessionLogs()->Display("*");


  //  chain->Process(pr,option.Data(),100,0);    
}

#ifndef __CINT__
//main function
int main(int argc, char** argv){
  theApp = new TRint("App",&argc,argv,NULL,0);
  //  theApp->Connect("keypressed(Int_t)","TSystem",gSystem,"ExitLoop()");
  //  TChain *t=new TChain("dstree");
  //  Readdatafile(t,Volume);
  // t->Add("/cache/shared/darkside/hqian/TMBfraction/outtmb_v8.root");
  int start, end;
  //read the input number from terminal
  if(theApp->Argc() == 2)
    {
      start = atoi(theApp->Argv(1));
      end = start;
    }
  else if(theApp->Argc() == 3)
    {
      start = atoi(theApp->Argv(1));
      end = atoi(theApp->Argv(2));
    }
  else{
    cout<<"Usage: ./reconMain startfile endfile "<<endl;
    cout<<"Usage: ./reconMain startfile "<<endl;
    return 0;
  }
  
  cout<<"Start Using the Selector..."<<endl;

  //read the dstree root file
  //#define use_outdir
#ifdef use_outdir
  string outdir = "/cache/shared/darkside/hqian/TMBfraction/cluster_data/";
#endif
  //  string output = outdir+ label+"_cylinder_clustered.root";
  for(int i=start; i<=end; i++)
    {
      string filename = Readdatafile(i);
      //   TFile* file = new TFile(filename.c_str());
      //   TTree* t = (TTree*) file->Get("dstree");
      TChain *t = new TChain("dstree");
      t->Add(filename.c_str());
#ifdef test_data
      string label ="out60gamma";
      //      string label = "outneutron";
#else
      string label = "outtmb";
#endif
      stringstream oss;
      if(i==0)
	label += "";
      else
	{
	  oss<<i;
	  label +="_v"+oss.str();
	}
#ifdef use_outdir
      label = outdir+label+"_clustered.root";
#else
      label +=  "_clustered.root";
#endif
      Process(t,label);
      oss.str("");
    }

  /*
   TChain *t = new TChain("dstree");
   string label = "outtmb";
   stringstream oss;
   Readdatafile(t,start,end);
   oss<<start<<end;
   label +="_v"+oss.str();
   label +=  "_cylinder_clustered.root";
   Process(t,label);
   oss.str("");
  */
 
 cout<<"Successfully finished the app."<<endl;
 return 1;

}
#endif /*__CINT__*/
