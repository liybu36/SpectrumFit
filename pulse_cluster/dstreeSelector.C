#define dstreeSelector_cxx
// The class definition in dstreeSelector.h has been generated automatically
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
// Root > T->Process("dstreeSelector.C")
// Root > T->Process("dstreeSelector.C","some options")
// Root > T->Process("dstreeSelector.C+")
//

#include "dstreeSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <vector>

using namespace std;
void dstreeSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void dstreeSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   SetOutputTree();

}

Bool_t dstreeSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either dstreeSelector::GetEntry() or TBranch::GetEntry()
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
  Int_t chainentry = fChain->GetChainEntryNumber(entry);
  fChain->GetEntry(entry);
  FillHistograms(entry);
  return kTRUE;
}

void dstreeSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void dstreeSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  
  string label = GetOption();
  string outdir = "/cache/shared/darkside/hqian/TMBfraction/cluster_data/";
  string output = outdir+ label+"_cylinder_clustered.root";

  TList* list = GetOutputList();
  ReconTree = dynamic_cast<TTree*>(list->FindObject(Form("ReconTree")));

  //Write histograms into a root file                                                                                               
  TFile f(output.c_str(), "RECREATE");
  ReconTree->Write();
  f.Write();
  f.Close();

}

void dstreeSelector::SetOutputTree()
{
  ReconTree = new TTree("Recon","Reconstructed Events");
  ReconTree->Branch("parent_energy", &recon_parent_enrg,32000,1);
  ReconTree->Branch("parent_pdg", &recon_parent_pdg,32000,1);
  ReconTree->Branch("event_id",&event_id);
  ReconTree->Branch("volume", &recon_vol, 32000, 1);
  ReconTree->Branch("particle", &recon_particle, 32000, 1);
  ReconTree->Branch("n_active",&n_active);
  ReconTree->Branch("et", &recon_et, 32000, 1);
  ReconTree->Branch("ex", &recon_ex, 32000, 1);
  ReconTree->Branch("ey", &recon_ey, 32000, 1);
  ReconTree->Branch("ez", &recon_ez, 32000, 1);
  ReconTree->Branch("dt", &recon_dt, 32000, 1);
  ReconTree->Branch("dz", &recon_dz, 32000, 1);
  ReconTree->Branch("dr", &recon_dr, 32000, 1);
  ReconTree->Branch("edep", &recon_edep, 32000, 1);
  ReconTree->Branch("edep_nuclear", &recon_edep_nuclear, 32000, 1);
  ReconTree->Branch("edep_electron", &recon_edep_electron, 32000, 1);
  ReconTree->Branch("eqch", &recon_eqch, 32000, 1);
  ReconTree->Branch("quenchingfactor", &recon_quenchingfactor, 32000, 1);
  //  ReconTree->Branch("edep_alpha", &recon_edep_alpha, 32000, 1);
  //  ReconTree->Branch("eqch_alpha", &recon_eqch_alpha, 32000, 1);
  ReconTree->Branch("contain_alpha", &recon_contain_alpha, 32000, 1);
  
  GetOutputList()->Add(ReconTree);
}

void dstreeSelector::InitializeRE(ReconEvent &evt, G4DSHitData &hit, int currentNum)
{
  evt.parent_en = hit.parent_enrg;
  evt.parent_pdg = hit.parent_pdg;
  evt.et = hit.time.at(currentNum);
  evt.ex = hit.x.at(currentNum);
  evt.ey = hit.y.at(currentNum);
  evt.ez = hit.z.at(currentNum);
  if(hit.volume.at(currentNum) =="p_scint")
    {
      evt.dt = veto_dt_cluster;
      evt.dz = veto_dz_cluster;
      evt.dr = veto_dr_cluster;
    }
  else
    {
      evt.dt = argon_dt_cluster;
      evt.dz = argon_dz_cluster;
      evt.dr = argon_dr_cluster;
      //      std::cout << "Z = " << evt.ez << endl;                                                                                
    }

  evt.particle = hit.particle.at(currentNum);
  evt.type = hit.type.at(currentNum);
  evt.edep = hit.enrg.at(currentNum);

  if (hit.particle.at(currentNum) == "alpha")
    {
      evt.edep_alpha = hit.enrg.at(currentNum);
      evt.contain_alpha = true;
      //      evt.eqch_alpha = hit.enrg.at(currentNum);                                                                             
    }
  else evt.contain_alpha = false;

  evt.eqch = hit.enrg.at(currentNum);//Initially set quenching to 1                                                                 
  evt.volume = hit.volume.at(currentNum);
  //  return;
}

bool dstreeSelector::CompatibleRE(ReconEvent evt, ReconEvent hit)
{

  if(hit.volume!= evt.volume || fabs(hit.et-evt.et) > evt.dt || fabs(hit.ez-evt.ez)>evt.dz)
    return false;
  if((hit.ex-evt.ex)*(hit.ex-evt.ex)+(hit.ey-evt.ey)*(hit.ey-evt.ey)>evt.dr*evt.dr)
    return false;

  return true;
}

void dstreeSelector::MergeRE(ReconEvent &evt, ReconEvent hit)
{
  Double_t dist, temp;
  evt.edep += hit.edep;
  evt.eqch += hit.eqch;

  if (hit.particle == "lepton" || hit.particle == "gamma")
    evt.edep_electron += hit.edep;
  if (hit.particle == "proton" || hit.particle == "neutron" || hit.particle == "alpha" || hit.particle == "nucleus")
    evt.edep_nuclear += hit.edep;
  if (hit.particle == "alpha")
    evt.contain_alpha = true;
  //  evt.edep_alpha += hit.edep_alpha;
  //  evt.eqch_alpha += hit.eqch_alpha;
  //  Double_t ratio = hit.edep/evt.edep;//fraction of energy increasement                                                          
  Double_t ratio = hit.eqch/evt.eqch;//should use quenched energy                                                                   

  //calculate the new barycenter of the time                                                                                        
  if(false && evt.volume == "p_scint"){
    if((hit.eqch > 2 && hit.et < evt.et) || evt.eqch < 2)
      evt.et = hit.et;
  }
  else{
    Double_t dist = fabs(evt.et-hit.et);
    evt.et += (hit.et-evt.et)*ratio;
    Double_t temp=dist*ratio;//distance to evt                                                                                      
    evt.dt += temp;
    temp = hit.dt+dist-temp;
    if(evt.dt<temp)evt.dt=temp;
  }
  //calculate the new barycenter of the z position                                                                                  
  dist = fabs(evt.ez-hit.ez);
  evt.ez += (hit.ez-evt.ez)*ratio;
  temp = dist*ratio;//distance to evt                                                                                               
  evt.dz += temp;
  temp = hit.dz+dist-temp;
  if(evt.dz<temp)
    evt.dz=temp;

  //calculate the new barycenter in the xy plane                                                                                    
  dist = sqrt((hit.ex-evt.ex)*(hit.ex-evt.ex)+(hit.ey-evt.ey)*(hit.ey-evt.ey));
  evt.ex += (hit.ex-evt.ex)*ratio;
  evt.ey += (hit.ey-evt.ey)*ratio;
  temp = dist*ratio;//distance to evt                                                                                               
  evt.dr += temp;
  temp = hit.dr+dist-temp;
  if(evt.dr<temp)
    evt.dr=temp;

  if(evt.particle!=hit.particle)
    evt.particle="mixed";
  if(evt.type!=hit.type)
    evt.type="mixed";

  return;
}

double dstreeSelector::get_point(vector<double> xs, vector<double> ys, double x)
{
  double deltay;
  double deltax;
  if(x < xs.at(0))
    {
      //cout << "WARNING : X = " << x                                                                                               
      //   << " BELOW RANGE, USING X(0) = " << xs.at(0) << endl;                                                                    
      return xs.at(0);
    }
  if(x > xs.at(xs.size()-1))
    {
      //cout << "WARNING : X = " << x                                                                                               
      //<< " ABOVE RANGE, USING X(" << xs.size()-1 <<") = " << xs.at(xs.size()-1)                                                   
      //<< endl;                                                                                                                    
      return xs.at(xs.size()-1);
    }
  for(unsigned int i = 0; i < xs.size(); i++)
    {
      if(xs.at(i) == x)
	return ys.at(i);
      else if(xs.at(i) > x)
        {
	  deltay = ys.at(i)-ys.at(i-1);
	  deltax = xs.at(i)-xs.at(i-1);
	  return ys.at(i-1)+(deltay/deltax)*(x-xs.at(i-1));
        }
    }
  // If code gets to here, something went wrong                                                                                     
  cout << "WARNiNG : SOMETHING WENT WRONG IN INTERPOLATION" << endl;
  return NULL;
}

double dstreeSelector::get_point(double *xs, double *ys, double x, int size)
{
  double deltay;
  double deltax;
  if(x < xs[0])
    {
      //      cout << "WARNING : X = " << x << " BELOW RANGE, USING X(0) = " << xs[0] << endl;                                      
      return ys[0];
    }
  if(x > xs[size-1])
    {
      //      cout << "WARNING : X = " << x << " ABOVE RANGE, USING X(" << size-1 <<") = " << xs[size-1] << endl;                   
      return ys[size-1];
    }
  for(int i = 0; i < size; i++)
    {
      if(xs[i] == x)
	return ys[i];
      else if(xs[i] > x)
        {
	  deltay = ys[i]-ys[i-1];
	  deltax = xs[i]-xs[i-1];
	  return ys[i-1]+(deltay/deltax)*(x-xs[i-1]);
        }
    }
  // If code gets to here, something went wrong                                                                                     
  cout << "WARNiNG : SOMETHING WENT WRONG IN INTERPOLATION" << endl;
  return NULL;
}

const int _ds10_quenching_ee_size = 251;
double _ds10_quenching_ee_energy_array[_ds10_quenching_ee_size] = //keVee                                                          
  {
    3, 5, 7, 9, 11,
    13, 15, 17, 19, 21,
    23, 25, 27, 29, 31,
    33, 35, 37, 39, 41,
    43, 45, 47, 49, 51,
    53, 55, 57, 59, 61,
    63, 65, 67, 69, 71,
    73, 75, 77, 79, 81,
    83, 85, 87, 89, 91,
    93, 95, 97, 99, 101,
    103, 105, 107, 109, 111,
    113, 115, 117, 119, 121,
    123, 125, 127, 129, 131,
    133, 135, 137, 139, 141,
    143, 145, 147, 149, 151,
    153, 155, 157, 159, 161,
    163, 165, 167, 169, 171,
    173, 175, 177, 179, 181,
    183, 185, 187, 189, 191,
    193, 195, 197, 199, 201,
    203, 205, 207, 209, 211,
    213, 215, 217, 219, 221,
    223, 225, 227, 229, 231,
    233, 235, 237, 239, 241,
    243, 245, 247, 249, 251,
    253, 255, 257, 259, 261,
    263, 265, 267, 269, 271,
    273, 275, 277, 279, 281,
    283, 285, 287, 289, 291,
    293, 295, 297, 299, 301,
    303, 305, 307, 309, 311,
    313, 315, 317, 319, 321,
    323, 325, 327, 329, 331,
    333, 335, 337, 339, 341,
    343, 345, 347, 349, 351,
    353, 355, 357, 359, 361,
    363, 365, 367, 369, 371,
    373, 375, 377, 379, 381,
    383, 385, 387, 389, 391,
    393, 395, 397, 399, 401,
    403, 405, 407, 409, 411,
    413, 415, 417, 419, 421,
    423, 425, 427, 429, 431,
    433, 435, 437, 439, 441,
    443, 445, 447, 449, 451,
    453, 455, 457, 459, 461,
    463, 465, 467, 469, 471,
    473, 475, 477, 479, 481,
    483, 485, 487, 489, 491,
    493, 495, 497, 499, 501,
    503
  };

double _ds10_quenching_ee_array[_ds10_quenching_ee_size] =
  {
    // Quenching at 1 kV/cm                                                                                                          
    // relative to light yield at null field (7.00 pe/keV)                                                                           
    // Values calculated from Run 3387                                                                                               
    0.596837, 0.643756, 0.664081, 0.679532, 0.681831,
    0.684578, 0.677978, 0.676092, 0.667289, 0.665633,
    0.651028, 0.648277, 0.648185, 0.640693, 0.635439,
    0.626229, 0.62093, 0.616566, 0.612934, 0.604095,
    0.604451, 0.6002, 0.603584, 0.591602, 0.586994,
    0.585982, 0.583503, 0.574758, 0.585129, 0.578986,
    0.569, 0.571115, 0.561534, 0.559565, 0.558337,
    0.557009, 0.55346, 0.556781, 0.551196, 0.550048,
    0.549547, 0.546511, 0.542514, 0.540375, 0.543217,
    0.540813, 0.533782, 0.535696, 0.536314, 0.532189,
    0.525858, 0.53233, 0.533782, 0.529218, 0.530761,
    0.524516, 0.53113, 0.526669, 0.522607, 0.514391,
    0.517642, 0.523029, 0.515696, 0.510605, 0.518794,
    0.511742, 0.510632, 0.515085, 0.517696, 0.512989,
    0.509259, 0.506989, 0.505043, 0.506972, 0.506673,
    0.504647, 0.50975, 0.503773, 0.503683, 0.501111,
    0.503756, 0.497568, 0.503536, 0.502822, 0.500786,
    0.504644, 0.492593, 0.50074, 0.497586, 0.494647,
    0.497549, 0.498107, 0.496273, 0.487788, 0.494165,
    0.490703, 0.486941, 0.484162, 0.490664, 0.49306,
    0.483707, 0.490684, 0.488487, 0.488284, 0.484356,
    0.487851, 0.487189, 0.479874, 0.483399, 0.483211,
    0.482808, 0.480194, 0.483197, 0.48607, 0.487002,
    0.480151, 0.483436, 0.481228, 0.481303, 0.479097,
    0.473938, 0.485983, 0.481064, 0.479787, 0.476905,
    0.480984, 0.478078, 0.475637, 0.476829, 0.478707,
    0.476295, 0.472569, 0.473348, 0.47567, 0.471471,
    0.477432, 0.475557, 0.474041, 0.473504, 0.478575,
    0.468693, 0.471295, 0.473281, 0.475602, 0.47028,
    0.471047, 0.47493, 0.468002, 0.472525, 0.473026,
    0.471801, 0.467276, 0.464736, 0.470248, 0.467677,
    0.465051, 0.473718, 0.468046, 0.469898, 0.470221,
    0.463342, 0.470322, 0.468003, 0.467527, 0.466834,
    0.466997, 0.464236, 0.465689, 0.467094, 0.463658,
    0.465637, 0.459848, 0.469017, 0.460924, 0.462493,
    0.469214, 0.466321, 0.466498, 0.467985, 0.459879,
    0.463663, 0.464527, 0.460028, 0.462595, 0.46456,
    0.464929, 0.471551, 0.461957, 0.464009, 0.461038,
    0.461585, 0.464325, 0.466728, 0.464421, 0.4653,
    0.462946, 0.465274, 0.459551, 0.462615, 0.466765,
    0.46087, 0.463251, 0.459242, 0.462542, 0.458255,
    0.467074, 0.46039, 0.462303, 0.457162, 0.463729,
    0.453992, 0.457272, 0.460397, 0.465321, 0.458937,
    0.453571, 0.460249, 0.457358, 0.456235, 0.464488,
    0.460977, 0.455952, 0.458472, 0.459481, 0.46219,
    0.458884, 0.460589, 0.464438, 0.455379, 0.459547,
    0.458208, 0.450635, 0.462524, 0.45625, 0.45676,
    0.460786, 0.465022, 0.455907, 0.460115, 0.460758,
    0.453698, 0.452668, 0.459375, 0.468595, 0.45782,
    0.473344, 0.453449, 0.456299, 0.46191, 0.464878,
    0.464878//Added so that extrapolation stays constant                                                                             
  };

//quenching data for protons in organic scintillator from Ben                                                                       
//AJW modified low energy points to agree with Astropart. Phys. 16 (2002)                                                            
//333-338                                                                                                                            
double MeV = 1000.0;
double dstreeSelector::o_scint_e_quenching (double e/*keV*/)
{
  //RNS                                                                                                                              
  //kB = 0.012 cm/MeV from Borexino                                                                                                  
  //Birks Quenching parameterized as in                                                                                              
  //"The ionization quench factor in liquid-scintillation counting standardizations                                                  
  //Malonda, Carles                                                                                                                  
  //Applied Radiation and Isotopes 51 (1999) 183-188                                                                                 

  double A1 = 0.32903;
  double A2 = 0.14404;
  double A3 = 0.08059;
  double A4 = -0.04536E-3;
  double A5 = 0.15623;
  double A6 = 0.06611;
  return ( (A1 + A2*TMath::Log(e) + A3*TMath::Log(e)*TMath::Log(e) + A4*TMath::Log(e)*TMath::Log(e)*TMath::Log(e))/
           (1 + A5*TMath::Log(e) + A6*TMath::Log(e)*TMath::Log(e) + A4*TMath::Log(e)*TMath::Log(e)*TMath::Log(e))  );
}

double o_scint_p_quenching_energy[] = {0,0.029,0.094,0.2,0.34,
                                       0.52, 0.72, 0.94,2,3,
                                       4,6,10,20,30,
                                       40,60,100};
double o_scint_p_quenching[] = {0,0.003,0.005,0.009,0.020,
                                0.041,0.071,0.127,0.6,1,
                                1.6,3,6,13,20,
                                30,45,70};

//quenching for carbon nuclear recoils                                                                                               
//AW from Astropart. Phys. 16 (2002) 333-338 and NIM 33 (1965) 131-135                                                              
 
double o_scint_c_quenching_energy[] = {0,0.046, 0.111, 0.229,
                                       0.368, 0.500, 1.2, 2.0,
                                       3.0, 4.0, 5.0};
double o_scint_c_quenching[] = {0,0.0022, 0.0026, 0.0032,
                                0.0044, 0.005, 0.007, 0.011,
				0.018, 0.025, 0.035};

double dstreeSelector::get_quenching_factor_organic_scintillator(ReconEvent evt)
{
  string particle(evt.particle);
  string type(evt.type);
  string volume(evt.volume);
  double e_dep = evt.edep;
  double e_quenched;

  // Load quenching data for particles in the scintillator                                                                          
  vector <double> pscintE_ins;
  vector <double> pscintE_outs;
  vector <double> ascintE_ins;
  vector <double> ascintE_outs;
  vector <double> cscintE_ins;
  vector <double> cscintE_outs;

  pscintE_ins.assign(o_scint_p_quenching_energy,o_scint_p_quenching_energy+sizeof(o_scint_p_quenching_energy)/sizeof(double));
  pscintE_outs.assign(o_scint_p_quenching,o_scint_p_quenching+sizeof(o_scint_p_quenching)/sizeof(double));

  ascintE_ins.assign(o_scint_a_quenching_energy, o_scint_a_quenching_energy+sizeof(o_scint_a_quenching_energy)/sizeof(double));
  ascintE_outs.assign(o_scint_a_quenching, o_scint_a_quenching+sizeof(o_scint_a_quenching)/sizeof(double));

  cscintE_ins.assign(o_scint_c_quenching_energy, o_scint_c_quenching_energy+sizeof(o_scint_c_quenching)/sizeof(double));
  cscintE_outs.assign(o_scint_c_quenching, o_scint_c_quenching+sizeof(o_scint_c_quenching)/sizeof(double));

  if(particle == "lepton")
    {
      return o_scint_e_quenching(e_dep);
      //return 1; //Paolo                                                                                                           
    }
  if(particle == "gamma")
    {
      // For now set to the same as electron quenching                                                                              
      return o_scint_e_quenching(e_dep);
      //      return 1; //Paolo                                                                                                     
    }
  if(type == "nucleus")
    {
      if(particle == "alpha")
	{
	  e_quenched = get_point(ascintE_ins, ascintE_outs, e_dep);
	  return 0.548*e_quenched/e_dep;
	  //      std::cout << "alpha : " << e_dep << " -> " << e_dep*0.027 << "\t OR "<< e_quenched << std::endl;                  
	  //return 0.027; // Paolo                                                                                                  
	}

      else
        {
          // Quench like carbon                                                                                                      
          e_quenched = get_point(cscintE_ins, cscintE_outs, e_dep);
          if(e_dep == 0)
            return 0;
          else
            return e_quenched/e_dep;
          //      return 0; // Paolo                                                                                                 
        }
    }
  if(particle == "proton")
    {

      e_quenched = get_point(pscintE_ins, pscintE_outs, e_dep);
      if(e_dep == 0)
	{
	  return 0;
	}
      else
	return e_quenched/e_dep;
      //return 0.333; // Paolo                                                                                                      
    }
  if(particle == "neutron")
    {
      e_quenched = get_point(pscintE_ins, pscintE_outs, e_dep);
      if(e_dep == 0)
	{
	  return 0;
	}
      else
	return e_quenched/e_dep;
      //return 0; // Paolo                                                                                                          
    }

  return 1;
}

double dstreeSelector::get_quenching_factor_liquid_argon(ReconEvent evt, bool fieldOn)
{
  string particle(evt.particle);
  string type(evt.type);
  string volume(evt.volume);
  double e_dep = evt.edep;

  // Load quenching data for particles in the scintillator                                                                           
  vector <double> pscintE_ins;
  vector <double> pscintE_outs;
  vector <double> ascintE_ins;
  vector <double> ascintE_outs;
  vector <double> cscintE_ins;
  vector <double> cscintE_outs;

  double QFn_e = 0.29;
  double QFn_f = 0.92;

  if(particle == "lepton")
    {
      if(fieldOn)
        return get_point(_ds10_quenching_ee_energy_array,
                         _ds10_quenching_ee_array,
                         e_dep,
                         _ds10_quenching_ee_size);
      else
        return 1;

    }
  if(particle == "gamma")
    {
      if(fieldOn) //For now use same quenching as for electrons                                                                     
        return get_point(_ds10_quenching_ee_energy_array,
                         _ds10_quenching_ee_array,
                         e_dep,
                         _ds10_quenching_ee_size);
      else
        return 1;

    }

  if(type == "nucleus")
    {
      if(particle == "alpha")
        {
          if(fieldOn)
            return QFn_e*QFn_f;
          else
            return QFn_e;
        }
      else
        {
          if(fieldOn)
            return QFn_e*QFn_f;
          else
            return QFn_e;
        }
    }
  if(particle == "proton")
    {
      if(fieldOn)
        return QFn_e*QFn_f;
      else
        return QFn_e;

    }

  if(particle == "neutron")
    {
      if(fieldOn)
        return QFn_e*QFn_f;
      else
        return QFn_e;
    }
  //  cout << "WARNING : PARTICLE " << particle << " of type " << type << " NOT FOUND. NOT QUENCHING" << endl;                       
  return 1;
}

double dstreeSelector::get_quenching_factor (ReconEvent evt, bool fieldOn)
{
  if (evt.volume == "p_active")
    return get_quenching_factor_liquid_argon(evt, fieldOn);
  else if (evt.volume == "p_scint")
    return get_quenching_factor_organic_scintillator(evt);
  else
    return 1.;
}

std::string dstreeSelector::convertPDG(int pdgNum){
  if(pdgNum == 22) return "gamma";
  else if(pdgNum == 11 || pdgNum == 12 || pdgNum == 13 || pdgNum == 14 || pdgNum == 15 || pdgNum == 16 || pdgNum == 17 || pdgNum == 18 ) return "lepton";
  else if(pdgNum == 2112) return "neutron";
  else if(pdgNum == 2212) return "proton";
  else if(pdgNum == 1000020040) return "alpha";
  else if(pdgNum > 1000010000) return "nucleus";
  else
    return "other";
}

bool dstreeSelector::checkActive(Float_t r, Float_t z){
  if( (r < maxR) && (minZ < z) && (z < maxZ) ){
    return true; //hit is within teflon cage                                                                                         
  }
  return false;
}

std::string dstreeSelector::convertVolume(int mat , Float_t r, Float_t z){
  if(mat == 2) return "p_scint";
  if( (mat == 8) && checkActive(r, z) ){
    return "p_active";
  }

  return "who_cares";
}

std::string dstreeSelector::returnType(int pdgNum){
  if(pdgNum == 22) return "boson";
  else if(pdgNum == 11 || pdgNum == 12 || pdgNum == 13 || pdgNum == 14 || pdgNum == 15 || pdgNum == 16 || pdgNum == 17 || pdgNum == 18 ) return "lepton";
  else if(pdgNum == 2112) return "baryon";
  else if(pdgNum == 2212) return "baryon";
  else if(pdgNum > 1000010000) return "nucleus";

  return "other";
}

void dstreeSelector::EmptyHit(G4DSHitData &currentHit)
{
  currentHit.x.clear();
  currentHit.y.clear();
  currentHit.z.clear();
  currentHit.r.clear();
  currentHit.enrg.clear();
  currentHit.time.clear();
  currentHit.particle.clear();
  currentHit.volume.clear();
  currentHit.type.clear();
  currentHit.deposits = 0;
  currentHit.parent_enrg = 0;
}

void dstreeSelector::returnHit(G4DSHitData &currentHit)
{
  EmptyHit(&currentHit);
  currentHit.parent_enrg=ene0*1000;
  currentHit.deposits = ndeposits;
  for(int j = 0; j < ndeposits; j++)
    {
      currentHit.x.push_back(dep_x[j]);
      currentHit.y.push_back(dep_y[j]);
      currentHit.z.push_back(dep_z[j]);
      currentHit.r.push_back(sqrt(dep_x[j]*dep_x[j] + dep_y[j]*dep_y[j])); 
      currentHit.particle.push_back(convertPDG( std::abs(dep_pdg[j]) ));
      currentHit.type.push_back(returnType( std::abs(dep_pdg[j]) ));
      currentHit.volume.push_back(convertVolume( dep_mat[j] , currentHit.r.at(j) , currentHit.z.at(j) ));
      currentHit.enrg.push_back(dep_ene[j] *1000); //convert meV to keV                
      currentHit.time.push_back(dep_time[j] * 1e-9 ); //convert ns to seconds                      
    }
}

void dstreeSelector::FillHistograms(Long64_t event)
{
  //  currentHit.returnHit(Events, event);  
  returnHit(&currentHit);
  n_active = 0;
  event_id = event;

  //clear all vectors for new event                                                                                               
  if(recon_evts->size())
    {
      recon_parent_enrg->clear();
      recon_parent_pdg->clear();
      recon_evts->clear();
      recon_et->clear();
      recon_ex->clear();
      recon_ey->clear();
      recon_ez->clear();
      recon_dt->clear();
      recon_dz->clear();
      recon_dr->clear();
      recon_edep->clear();
      recon_edep_electron->clear();
      recon_edep_nuclear->clear();
      recon_eqch->clear();
      recon_vol->clear();
      recon_particle->clear();
      recon_quenchingfactor->clear();
      //      recon_edep_alpha->clear();
      //      recon_eqch_alpha->clear();
      recon_contain_alpha->clear();
    }

  for(size_t j=0; j<(size_t)currentHit.deposits; j++)
    {

      if(currentHit.volume.at(j)=="who_cares" || currentHit.enrg.at(j)<1e-3)
	continue;//remove 0 energy events                                                                                         

      ReconEvent anEvt;
      InitializeRE(anEvt,currentHit , j);

      double quenchfactor = get_quenching_factor(anEvt, true);
      anEvt.eqch = anEvt.edep * quenchfactor;

      //      if(currentHit.particle.at(j)=="alpha")                                                                              
      //      if(anEvt.particle == "alpha")
      //	anEvt.eqch_alpha = anEvt.edep_alpha * quenchfactor;

      bool chained = false;

      //loop over existing clusters in the current event                                                                          
      for(size_t k=0; k<recon_evts->size(); k++)
	{
	  //Check if current hit should be included in an existing cluster                                                        
	  if(CompatibleRE(recon_evts->at(k),anEvt))
	    {
	      MergeRE(recon_evts->at(k),anEvt);
	      chained = true;
	      break;
	    }
	}

      if(!chained)
	recon_evts->push_back(anEvt);

    }//end for j, scan all hits                                                                                                   

  //if no active energy deposition                                                                                              
  if(recon_evts->size()<1)
    continue;

  //need to delete events with volume=-1                                                                                        
  for(size_t k=0; k<recon_evts->size(); k++)
    {
      if(recon_evts->at(k).volume=="who_cares")
	{
	  recon_evts->erase(recon_evts->begin()+k);
	  k--;
	}
    }

  //need to sort the events by time                                                                                             
  if(recon_evts->size()>1)
    std::sort(recon_evts->begin(), recon_evts->end(), CompareRE);
  //Fill output TTree variables                                                                                                 
  for(size_t k=0; k<recon_evts->size(); k++)
    {
      if(recon_evts->at(k).volume == "who_cares")
	continue; //has been deleted                                                                                            
      //need energy cut 3keVee?                                                                                                 
      if(recon_evts->at(k).volume=="p_active")
	{  n_active++;
	}
      recon_evts->at(k).quenchingfactor = recon_evts->at(k).eqch/recon_evts->at(k).edep;
      //      std::cout << recon_evts->at(k).parent_en << std::endl;                                                            
      recon_parent_enrg->push_back(recon_evts->at(k).parent_en);
      recon_parent_pdg->push_back(recon_evts->at(k).parent_pdg);
      recon_et->push_back(recon_evts->at(k).et);
      recon_ex->push_back(recon_evts->at(k).ex);
      recon_ey->push_back(recon_evts->at(k).ey);
      recon_ez->push_back(recon_evts->at(k).ez);
      recon_dt->push_back(recon_evts->at(k).dt);
      recon_dz->push_back(recon_evts->at(k).dz);
      recon_dr->push_back(recon_evts->at(k).dr);
      recon_edep->push_back(recon_evts->at(k).edep);
      recon_edep_nuclear->push_back(recon_evts->at(k).edep_nuclear);
      recon_edep_electron->push_back(recon_evts->at(k).edep_electron);
      recon_eqch->push_back(recon_evts->at(k).eqch);
      recon_quenchingfactor->push_back(recon_evts->at(k).quenchingfactor);
      recon_vol->push_back(recon_evts->at(k).volume);
      recon_particle->push_back(recon_evts->at(k).particle);
      //    recon_edep_alpha->push_back(recon_evts->at(k).edep_alpha);
      //    recon_eqch_alpha->push_back(recon_evts->at(k).eqch_alpha);
      recon_contain_alpha->push_back(recon_evts->at(k).contain_alpha);
    }//end k                                      
  if(recon_evts->size())
    ReconTree->Fill();
}
