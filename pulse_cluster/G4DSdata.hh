#ifndef G4DSDATA_hh
#define G4DSDATA_hh

#include <TObject.h>
#include "TROOT.h"
#include "TTree.h"
#include "TObject.h"
#include <string>
#include <vector>
#include <cstdlib>

typedef std::vector<Double_t>  dvec;
typedef std::vector<Float_t> fvec;
typedef std::vector< std::string > svec;

//struct G4DSHitData : public TObject
struct G4DSHitData
{
  G4DSHitData() :
    deposits(0), parent_pdg(0),
    //    x.clear(),y.clear(),z.clear(),r.clear(),enrg.clear(),time.clear(),
    //   particle.clear(), volume.clear(), type.clear()
    x(0),y(0),z(0),r(0),parent_enrg(0),
    enrg(0),time(0),particle(0),volume(0),type(0)
  {}
  //  virtual ~G4DSHitData();

  int deposits;
  int parent_pdg;
 
  fvec x;
  fvec y;
  fvec z;
  fvec r;

  Float_t parent_enrg;  
  fvec enrg;
  dvec time;
  
  svec particle;
  svec volume;
  svec type;

  //  void Empty();
  //  void resize(int size);
  //  void returnHit(TTree *tree, int entryNum);

  //  ClassDef(G4DSHitData, 1)
};

struct ReconEvent
{
  ReconEvent():
    parent_pdg(0), et(0), ex(0), ey(0), ez(0), dt(0), dz(0), dr(0),edep(0),edep_nuclear(0),
    edep_electron(0),eqch(0),quenchingfactor(0),volume(0), particle(0),type(0),parent_en(0),
    contain_alpha(0)
  {}
  
  //need another member to record fraction of quenching                                                                             
  Int_t parent_pdg;
  Double_t et;
  Double_t ex;
  Double_t ey;
  Double_t ez;
  Double_t dt;
  Double_t dz;
  Double_t dr;
  Double_t edep;
  Double_t edep_nuclear;
  Double_t edep_electron;
  Double_t eqch;
  Double_t quenchingfactor;
  TString volume;
  TString particle;
  TString type;//not saved in root tree                                                                                              
  Float_t parent_en;

  //  Double_t edep_alpha;
  //  Double_t eqch_alpha;
  Bool_t contain_alpha;

  void clear(){
    parent_pdg=0; et=0; ex=0; ey=0; ez=0; dt=0; dz=0; dr=0; edep=0; edep_nuclear=0;
    edep_electron=0; eqch=0; quenchingfactor=0; parent_en=0;
    contain_alpha=false;  
    volume=""; particle=""; type="";
  }

};


#endif
