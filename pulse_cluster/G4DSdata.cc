#include "G4DSdata.hh"
#include <iostream>
#include <stdio.h>
#include <math.h>

using namespace std;
int MAXDEP = 70000; //from g4rooter
Float_t zshift = -1;
Float_t maxR = 17.77;  //all in cm
Float_t maxZ = 14.6; 
Float_t minZ = -22.2;
//Float_t maxZ = 27.5; //=center+halfz of teflon_h  
//Float_t minZ = -8.6;   //=center-halfz of teflon_h
//center=91.81cm
//halfz of teflon_h= G2: 67.5cm
//radius of DSG2=74.15cm
//halfz of DS50=19.304cm
//radius of DS50=17.77cm

G4DSHitData::G4DSHitData(){
}

G4DSHitData::~G4DSHitData(){
}

std::string convertPDG(int pdgNum){
  if(pdgNum == 22) return "gamma";
  else if(pdgNum == 11 || pdgNum == 12 || pdgNum == 13 || pdgNum == 14 || pdgNum == 15 || pdgNum == 16 || pdgNum == 17 || pdgNum == 18 ) return "lepton";
  else if(pdgNum == 2112) return "neutron";
  else if(pdgNum == 2212) return "proton";
  else if(pdgNum == 1000020040) return "alpha";
  else if(pdgNum > 1000010000) return "nucleus";
  else
  return "other"; 
}

bool checkActive(Float_t r, Float_t z){
  if( (r < maxR) && (minZ < z) && (z < maxZ) ){
    return true; //hit is within teflon cage
  }
  return false;
}

std::string convertVolume(int mat , Float_t r, Float_t z){
  if(mat == 2) return "p_scint";
  if( (mat == 8) && checkActive(r, z) ){
    return "p_active";
  }

  return "who_cares";
}

std::string returnType(int pdgNum){
  if(pdgNum == 22) return "boson";
  else if(pdgNum == 11 || pdgNum == 12 || pdgNum == 13 || pdgNum == 14 || pdgNum == 15 || pdgNum == 16 || pdgNum == 17 || pdgNum == 18 ) return "lepton";
  else if(pdgNum == 2112) return "baryon";
  else if(pdgNum == 2212) return "baryon";
  else if(pdgNum > 1000010000) return "nucleus";
  
  return "other"; 
}

void G4DSHitData::Empty()
{
  x.clear();
  y.clear();
  z.clear();
  r.clear();
  enrg.clear();
  time.clear();
  particle.clear();
  volume.clear();
  type.clear();  
  deposits=0;
  parent_enrg=0;
}

void G4DSHitData::resize(int size)
{
  x.resize(size);
  y.resize(size);
  z.resize(size);
  r.resize(size);
  enrg.resize(size);
  time.resize(size);
  particle.resize(size);
  volume.resize(size);  
  type.resize(size);
}

void G4DSHitData::returnHit(TTree *tree, int entryNum ){
  Empty();
  resize(1);
  int ndep;
  Float_t dep_e[MAXDEP], dep_x[MAXDEP], dep_y[MAXDEP], dep_z[MAXDEP], dep_r[MAXDEP];
  Int_t dep_pdg[MAXDEP], dep_mat[MAXDEP];
  Double_t dep_t[MAXDEP];
  Float_t ene0;

  tree->SetBranchAddress("ndeposits", &ndep);
  tree->SetBranchAddress("dep_ene", &dep_e);
  tree->SetBranchAddress("dep_x", &dep_x);
  tree->SetBranchAddress("dep_y", &dep_y);
  tree->SetBranchAddress("dep_z", &dep_z);
  tree->SetBranchAddress("dep_r", &dep_r);
  tree->SetBranchAddress("dep_pdg", &dep_pdg);
  tree->SetBranchAddress("dep_mat", &dep_mat);
  tree->SetBranchAddress("dep_time", &dep_t);
  tree->SetBranchAddress("ene0", &ene0);

  tree->GetEntry(entryNum);
  parent_enrg=ene0*1000;
  for(int j = 0; j < ndep; j++){
    resize(ndep);
    x.at(j) = dep_x[j];
    y.at(j) = dep_y[j];
    z.at(j) = dep_z[j];
    r.at(j) = sqrt(dep_x[j]*dep_x[j] + dep_y[j]*dep_y[j]);//dep_r[j];
    particle.at(j) = convertPDG( std::abs(dep_pdg[j]) );
    type.at(j) = returnType( std::abs(dep_pdg[j]) );
    //    if(dep_mat[j] == 8 && dep_z[j] > 15)
    //      std::cout << "Z = " << z.at(j) << " : " << dep_z[j] << "\t r = " << dep_r[j] << " : " << sqrt(dep_x[j]*dep_x[j]+dep_y[j]*dep_y[j]) << std::endl;
    volume.at(j) = convertVolume( dep_mat[j] , r.at(j) , z.at(j) );
    enrg.at(j) = dep_e[j] *1000; //convert meV to keV
    time.at(j) = dep_t[j] * 1e-9 ; //convert ns to seconds
    deposits = ndep;
  }
}
