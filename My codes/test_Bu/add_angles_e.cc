#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>
#include <cmath>
#include <TRandom.h>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <string>

using std::cout;
using TMath::Cos;
using vec4d = ROOT::Math::PxPyPzEVector;
using ROOT::Math::VectorUtil::boost;
using ROOT::Math::VectorUtil::Angle;
using vec3d =  ROOT::Math::XYZVector;
using namespace std;
void lhcb_angles(bool bzero, const vec4d & eplus, const vec4d & eminus, const vec4d & kaon, const vec4d & pion, double& costhetal, double& costhetak, double& phi);

template<typename T> struct BranchGetter
{
  T data ;
  operator T&() { return data ; }
  operator const T&() const { return data ; }
  BranchGetter(TTree& tree, const TString& fieldname, bool reportIfMissing = true) {
    tree.SetBranchAddress(fieldname,&data) ;
    tree.SetBranchStatus(fieldname,1) ;
  }
} ;

template<typename T> struct BranchSetter
{
  T data ;
  operator T&() { return data ; }
  BranchSetter(TTree& tree, const TString& fieldname) {
    tree.Branch(fieldname,&data) ;
  }
  BranchSetter& operator=(const T& x) { data = x ; return *this ;}
} ;


struct ParticleGetter
{
  BranchGetter<Int_t> ID ;
  BranchGetter<Double_t> PX ;
  BranchGetter<Double_t> PY ;
  BranchGetter<Double_t> PZ ;
  BranchGetter<Double_t> PE ;
  BranchGetter<Double_t> M ;
  ParticleGetter(TTree& tree, TString name)
    : ID{tree,name+"_ID"},
      PX{tree,name+"_PX"},
      PY{tree,name+"_PY"},
      PZ{tree,name+"_PZ"},
      PE{tree,name+"_PE"},
      M{tree,name+"_M"} {}
  auto p4() const { return ROOT::Math::PxPyPzEVector{PX,PY,PZ,PE} ; }
} ;


struct TrParticleGetter
{
  BranchGetter<Double_t> TrPX ;
  BranchGetter<Double_t> TrPY ;
  BranchGetter<Double_t> TrPZ ; 
  TrParticleGetter(TTree& tree, TString name)
    : TrPX{tree,name+"_TrPX"},
      TrPY{tree,name+"_TrPY"},
      TrPZ{tree,name+"_TrPZ"} {}
  auto p4(const double M) const { return ROOT::Math::PxPyPzMVector{TrPX,TrPY,TrPZ,M} ; }
} ;

struct ElectronParticleGetter : ParticleGetter
{
  BranchGetter<Bool_t>   HasBremAdded;
  BranchGetter<Double_t> BremMultiplicity; // why is this a double?
  BranchGetter<Double_t> BremPE;
  BranchGetter<Double_t> BremPX;
  BranchGetter<Double_t> BremPY;
  BranchGetter<Double_t> BremPZ;
  ElectronParticleGetter(TTree& tree, TString name)
    : ParticleGetter{tree,name},
      HasBremAdded{tree,name+"_HasBremAdded"},
      BremMultiplicity{tree,name+"_BremMultiplicity"},
      BremPE{tree,name+"_BremPE"},
      BremPX{tree,name+"_BremPX"},
      BremPY{tree,name+"_BremPY"},
      BremPZ{tree,name+"_BremPX"} {}
} ;

struct DTFParticleGetter
{
  //BranchGetter<Int_t> ID ;
  BranchGetter<Double_t> PX ;
  BranchGetter<Double_t> PY ;
  BranchGetter<Double_t> PZ ;
  Double_t mass2 ;
  DTFParticleGetter(TTree& tree, TString name, double m2)
    : PX{tree,name+"_PX_0"},
      PY{tree,name+"_PY_0"},
      PZ{tree,name+"_PZ_0"},
      mass2{m2}{}
  auto p4() const {
    auto PE = std::sqrt( mass2 + PX*PX + PY*PY + PZ*PZ ) ;
    return ROOT::Math::PxPyPzEVector{PX,PY,PZ,PE} ; }
} ;

struct TrueParticleGetter
{
  BranchGetter<Int_t> ID ;
  BranchGetter<Double_t> PX ;
  BranchGetter<Double_t> PY ;
  BranchGetter<Double_t> PZ ;
  BranchGetter<Double_t> PE ;
   TrueParticleGetter(TTree& tree, TString name)
    : ID{tree,name+"_TRUEID"},
      PX{tree,name+"_TRUEP_X"},
      PY{tree,name+"_TRUEP_Y"},
      PZ{tree,name+"_TRUEP_Z"},
      PE{tree,name+"_TRUEP_E"} {}
  auto p4() const { return ROOT::Math::PxPyPzEVector{PX,PY,PZ,PE} ; }
  //auto tp4() const { return TLorentzVector{PX,PY,PZ,PE} ; }
} ;

struct MCDT_TrueParticleGetter
{
  BranchGetter<Int_t> ID ;
  BranchGetter<Double_t> PX ;
  BranchGetter<Double_t> PY ;
  BranchGetter<Double_t> PZ ;
  BranchGetter<Double_t> PE ;
   MCDT_TrueParticleGetter(TTree& tree, TString name)
    : ID{tree,"MCDT_"+name+"_ID"},
      PX{tree,"MCDT_"+name+"_TRUEP_X"},
      PY{tree,"MCDT_"+name+"_TRUEP_Y"},
      PZ{tree,"MCDT_"+name+"_TRUEP_Z"},
      PE{tree,"MCDT_"+name+"_TRUEP_E"} {}
  auto p4() const { return ROOT::Math::PxPyPzEVector{PX,PY,PZ,PE} ; }
  //auto tp4() const { return TLorentzVector{PX,PY,PZ,PE} ; }
} ;

struct FourVectorSetter
{
  BranchSetter<Double_t> PX ;
  BranchSetter<Double_t> PY ;
  BranchSetter<Double_t> PZ ;
  BranchSetter<Double_t> PE ;
  FourVectorSetter(TTree& tree, const TString& name)
    : PX{tree, name + "P_X"},
      PY{tree, name + "P_Y"},
      PZ{tree, name + "P_Z"},
      PE{tree, name + "P_E"} {}
  template<typename T>
  FourVectorSetter& operator=( const T& p4 ) { PX = p4.x(); PY = p4.y() ; PZ = p4.z(); PE = p4.E() ; return *this ;}
  void reset() { PX=PY=PZ=PE=0 ; }
} ;



auto throwlandautruncated(double maxval=100.)
{
  // throw energyloss based on Landau, but truncate at maxval
  double rc{0} ;
  do {
    rc = gRandom->Landau(0,1) ;
  } while(rc>maxval) ;
  return rc ;
}

// Rewritten by WH, 10.11.21
void add_angles_e(TString const fileName, TString const treeName, TString const fileNameNew, bool use_MC, std::string const channel){
  
  /////////////////////////////////////////New M.Geijsen & A.Biolchini  Dec-jan 2022-2023
  const double mpion = 139.57 ;
  const double mkaon = 493.68 ;
  const double melec = 0.511 ;
  /////////////////////////////////////////END
      cout << "Processing MC" << endl;

  TString inputDir = "";
  TFile f(inputDir+fileName);
  if(! f.IsOpen()) {
    std::cout << "Cannot find input file: \" "
              << inputDir+fileName << "\"" << std::endl ;
    return ;
  }
  
  auto tree = dynamic_cast<TTree*>(f.Get(treeName)) ;
  if( !tree ) {
    std::cout << "Cannot find tree with name " << treeName << std::endl 
	      << "File contents are: " << std::endl;
    f.ls() ;
    return ;
  }
  
  std::cout << "Running add_angles_e.cc" << std::endl
	    << "File name (reading):" << std::endl
	    << fileName << std::endl
	    << "File name (saving to):" << std::endl
	    << fileNameNew << std::endl ;

  TFile *newfile = new TFile(inputDir+fileNameNew, "RECREATE");
  TTree *newtree = tree->CloneTree(0); // Do no copy the data yet



  //tree->SetBranchStatus("*",0) ;
  std::unique_ptr<ParticleGetter> kaon , pion, Kstar, B, JPsi;
  kaon.reset(new ParticleGetter {*tree,"K"} ) ;
  pion.reset(new ParticleGetter {*tree,"Pi"} );
  Kstar.reset(new ParticleGetter {*tree,"Kstar"} ) ;
  B.reset(new ParticleGetter {*tree,"B"} );
  JPsi.reset(new ParticleGetter {*tree,"Jpsi"} ) ;

  std::unique_ptr<ElectronParticleGetter> L1, L2;
  L1.reset(new ElectronParticleGetter {*tree,"L1"} );
  L2.reset(new ElectronParticleGetter {*tree,"L2"} );



  // These are missing from the RX ntuples so we cannot compute the angles!
  //constexpr Double_t me2 = 0.511 * 0.511;
  //constexpr Double_t mk2 = 493.6770041 * 493.6770041;
  //constexpr Double_t mpi2 = 139.57017584 * 139.57017584;
  //DTFParticleGetter dtfL1{*tree,"B_DTF_PVandB_L1",me2} ;
  //DTFParticleGetter dtfL2{*tree,"B_DTF_PVandB_L2",me2} ;
  //DTFParticleGetter dtfKaon{*tree,"B_DTF_PVandB_K",mk2} ;
  //DTFParticleGetter dtfPion{*tree,"B_DTF_PVandB_Pi",mpi2} ;
  //Double_t ctk_dtf, ctl_dtf, phi_dtf; 
  //newtree->Branch("ctk_DTF_PVandB_0", &ctk_dtf);
  //newtree->Branch("ctl_DTF_PVandB_0", &ctl_dtf);
  //newtree->Branch("phi_DTF_PVandB_0", &phi_dtf);
 
  Double_t ctk_lhcb, ctl_lhcb, phi_lhcb; 
  newtree->Branch("ctk_lhcb", &ctk_lhcb);
  newtree->Branch("ctl_lhcb", &ctl_lhcb);
  newtree->Branch("phi_lhcb", &phi_lhcb);

  /////////////////////////////////////////New M.Geijsen & A.Biolchini  Dec-jan 2022-2023
  //add new branches:
  Double_t B_M_SubstL1_e2pi, B_M_SubstL1L2_ee2pipi, B_M_SubstL1L2_ee2kk;

  newtree->Branch("B_M_SubstL1_e2pi", &B_M_SubstL1_e2pi);
  newtree->Branch("B_M_SubstL1L2_ee2pipi", &B_M_SubstL1L2_ee2pipi);
  newtree->Branch("B_M_SubstL1L2_ee2kk", &B_M_SubstL1L2_ee2kk);
 
  //Defining variable to fill only when sample is Data
  Double_t B_M_Tr_SubstL1_e2pi, B_M_Tr_SubstL1L2_ee2pipi, B_M_Tr_SubstL1L2_ee2kk;
  std::unique_ptr<TrParticleGetter> Trkaon, Trpion, TrL1, TrL2;
  /////////////////////////////////////////END
  
    // including recovered photons - M.Soares 20.Oct.21

    // A.Biolchini & M.Soares 3.Jan.2023
    // Commented out because the flag `reco_wrongID` has been added to keep trace of the B candidates that were not correctly reconstructed
    // All these true variable were not correct, because they were containing the information of particles not correctly reconstructed.
    //Double_t ctk_lhcb_true, ctl_lhcb_true, phi_lhcb_true;
    //Int_t nPhotonRecov ;
    //Double_t ctk_lhcb_recov, ctl_lhcb_recov, phi_lhcb_recov;
    //std::unique_ptr<FourVectorSetter>   L1_RECOVP, L2_RECOVP;

  std::unique_ptr<TrueParticleGetter> trueKaon,truePion,trueL1,trueL2,trueKstar,trueB ;
  Int_t reco_wrongID;

  Double_t MCDT_ctk_lhcb_true, MCDT_ctl_lhcb_true, MCDT_phi_lhcb_true;
  Int_t MCDT_nPhotonRecov ;
  Double_t MCDT_ctk_lhcb_recov, MCDT_ctl_lhcb_recov, MCDT_phi_lhcb_recov;
  std::unique_ptr<FourVectorSetter>   MCDT_L1_RECOVP, MCDT_L2_RECOVP;
  std::unique_ptr<MCDT_TrueParticleGetter> MCDT_trueKaon,MCDT_truePion,MCDT_trueL1,MCDT_trueL2,MCDT_trueKstar,MCDT_trueB ;

  std::unique_ptr<FourVectorSetter>   L1_SMEARP, L2_SMEARP;
  Double_t B_M_smeared, Jpsi_M_smeared;


  if(use_MC==true){
    cout << "Processing MC" << endl;
    
    // BranchSetter<Double_t> B_M_smeared(*newtree,"B_M_smeared") ;
    // BranchSetter<Double_t> Jpsi_M_smeared(*newtree,"Jpsi_M_smeared") ;
    newtree->Branch("B_M_smeared", &B_M_smeared);
    newtree->Branch("Jpsi_M_smeared", &Jpsi_M_smeared);
    L1_SMEARP.reset( new FourVectorSetter{*newtree, "L1_SMEAR"} ) ;
    L2_SMEARP.reset( new FourVectorSetter{*newtree, "L2_SMEAR"} ) ;

    newtree->Branch("reco_wrongID", &reco_wrongID);
    trueKaon.reset(new TrueParticleGetter{*tree,"K"}) ;
    truePion.reset(new TrueParticleGetter{*tree,"Pi"}) ;
    trueL1.reset(new TrueParticleGetter{*tree,"L1"}) ;
    trueL2.reset(new TrueParticleGetter{*tree,"L2"}) ;
    trueKstar.reset(new TrueParticleGetter{*tree,"Kstar"}) ;
    trueB.reset(new TrueParticleGetter{*tree,"B"}) ;

    //newtree->Branch("ctk_lhcb_true", &ctk_lhcb_true);
    //newtree->Branch("ctl_lhcb_true", &ctl_lhcb_true);
    //newtree->Branch("phi_lhcb_true", &phi_lhcb_true);

    // including recovered photons - M.Soares 20.Oct.21
    //newtree->Branch("ctk_lhcb_recov", &ctk_lhcb_recov);
    //newtree->Branch("ctl_lhcb_recov", &ctl_lhcb_recov);
    //newtree->Branch("phi_lhcb_recov", &phi_lhcb_recov);
    //newtree->Branch("nPhotonRecov", &nPhotonRecov);
    //L1_RECOVP.reset( new FourVectorSetter{*newtree, "L1_RECOV"} ) ;
    //L2_RECOVP.reset( new FourVectorSetter{*newtree, "L2_RECOV"} ) ;

    if (channel == "B02Kst0Jpsi2ee" || channel == "B02Kst0eeSignal")
    {
      MCDT_trueKaon.reset(new MCDT_TrueParticleGetter{*tree,"K"}) ;
      MCDT_truePion.reset(new MCDT_TrueParticleGetter{*tree,"Pi"}) ;
      MCDT_trueL1.reset(new MCDT_TrueParticleGetter{*tree,"L1"}) ;
      MCDT_trueL2.reset(new MCDT_TrueParticleGetter{*tree,"L2"}) ;
      MCDT_trueKstar.reset(new MCDT_TrueParticleGetter{*tree,"Kstar"}) ;
      MCDT_trueB.reset(new MCDT_TrueParticleGetter{*tree,"B"}) ;

      newtree->Branch("MCDT_ctk_lhcb_true", &MCDT_ctk_lhcb_true);
      newtree->Branch("MCDT_ctl_lhcb_true", &MCDT_ctl_lhcb_true);
      newtree->Branch("MCDT_phi_lhcb_true", &MCDT_phi_lhcb_true);

      newtree->Branch("MCDT_ctk_lhcb_recov", &MCDT_ctk_lhcb_recov);
      newtree->Branch("MCDT_ctl_lhcb_recov", &MCDT_ctl_lhcb_recov);
      newtree->Branch("MCDT_phi_lhcb_recov", &MCDT_phi_lhcb_recov);
      newtree->Branch("MCDT_nPhotonRecov", &MCDT_nPhotonRecov);
      MCDT_L1_RECOVP.reset( new FourVectorSetter{*newtree, "MCDT_L1_RECOV"} ) ;
      MCDT_L2_RECOVP.reset( new FourVectorSetter{*newtree, "MCDT_L2_RECOV"} ) ;

      std::cout << "MCDT angles only added for channels B02Kst0Jpsi2ee and B02Kst0eeSignal" << std::endl;
    }
  }
  else {
  
    std::cout << "Processing data" << std::endl;

    /////////////////////////////////////////New M.Geijsen & A.Biolchini  Dec-jan 2022-2023
    Trkaon.reset(new TrParticleGetter{*tree,"K"}) ;
    Trpion.reset(new TrParticleGetter{*tree,"Pi"}) ;
    TrL1.reset(new TrParticleGetter{*tree,"L1"}) ;
    TrL2.reset(new TrParticleGetter{*tree,"L2"}) ;

    // add new branches:
    newtree->Branch("B_M_Tr_SubstL1_e2pi", &B_M_Tr_SubstL1_e2pi);
    newtree->Branch("B_M_Tr_SubstL1L2_ee2pipi", &B_M_Tr_SubstL1L2_ee2pipi);
    newtree->Branch("B_M_Tr_SubstL1L2_ee2kk", &B_M_Tr_SubstL1L2_ee2kk);
    ////////////////////////////////////////END
  }

  int const Nevent=tree->GetEntries();
  int bin;

   // Define the true ids for each channel
   std::unordered_map<std::string, std::array<int, 5>> true_ids;
   true_ids["B02Kst0eeFlat"]    = {511, 321, 211, 11, 11};
   true_ids["B02Kst0eeSignal"]  = {511, 321, 211, 11, 11};
   true_ids["B02Kst0Jpsi2ee"]   = {511, 321, 211, 11, 11};
   true_ids["B02Kst0Psi2S2ee"]  = {511, 321, 211, 11, 11};
   true_ids["Bd2Denu2Kst0enu"]  = {511, 321, 211, 11, 11};
   true_ids["Bu2Kpipiee"]       = {521, 321, 211, 11, 11};
   true_ids["Lb02pKJpsi2ee"]    = {5122, 321, 2212, 11, 11};
   true_ids["Bu2KpipiPsi2S2ee"] = {521, 321, 211, 11, 11};


  for (int kk=0;kk<Nevent; ++kk){
    if(kk%10000==0) std:: cout<<kk<<" out of "<<Nevent<< std::endl;
    tree->GetEntry(kk);
    
    const auto eplus  = L1->ID<0 ? L1->p4() : L2->p4() ;
    const auto eminus = L1->ID<0 ? L2->p4() : L1->p4() ;
    auto bzero = kaon->ID > 0 ;
    
    /////////////////////////////////////////New M.Geijsen & A.Biolchini & M.Soares Dec-jan 2022-2023
    //
    ROOT::Math::PxPyPzMVector kaon_M, pion_M, L2_M ;

    kaon_M.SetCoordinates(kaon->p4().Px(),kaon->p4().Py(),kaon->p4().Pz(),mkaon) ;
    pion_M.SetCoordinates(pion->p4().Px(),pion->p4().Py(),pion->p4().Pz(),mpion) ;

    // We could have it here for completeness, but never used
    // ROOT::Math::PxPyPzMVector L1_M
    // L1_M.SetCoordinates(L1->p4().Px(),L1->p4().Py(),L1->p4().Pz(),melec) ;
    L2_M.SetCoordinates(L2->p4().Px(),L2->p4().Py(),L2->p4().Pz(),melec) ;

    // Changing masses to pion or kaon  
    ROOT::Math::PxPyPzMVector  L1_Mpi, L2_Mpi, L1_Mka, L2_Mka;
    L1_Mpi.SetCoordinates(L1->p4().Px(),L1->p4().Py(),L1->p4().Pz(),mpion) ;
    L2_Mpi.SetCoordinates(L2->p4().Px(),L2->p4().Py(),L2->p4().Pz(),mpion) ;

    L1_Mka.SetCoordinates(L1->p4().Px(),L1->p4().Py(),L1->p4().Pz(),mkaon) ;
    L2_Mka.SetCoordinates(L2->p4().Px(),L2->p4().Py(),L2->p4().Pz(),mkaon) ;



    B_M_SubstL1_e2pi = ( kaon_M +pion_M + L1_Mpi + L2_M ).M() ;
    B_M_SubstL1L2_ee2pipi = ( kaon_M +pion_M + L1_Mpi + L2_Mpi ).M() ;
    B_M_SubstL1L2_ee2kk = ( kaon_M +pion_M + L1_Mka +  L2_Mka ).M() ; 

    //cout<<"Check masses:"<< kaon_M.M() <<"  "<<kaon->p4().M() <<"  " << kaon->M << std::endl ;
    //
    ///////////////////////////////////////END


    lhcb_angles(bzero, eplus, eminus, kaon->p4(), pion->p4(), ctl_lhcb, ctk_lhcb, phi_lhcb);

    //debugging
    ////cout<<"Before :"<<ctl_lhcb<< std::endl;


    // Using DTF PEXYZ
    //const auto dtf_eplus  = L1.ID<0 ? dtfL1.p4() : dtfL2.p4() ;
    //const auto dtf_eminus = L1.ID<0 ? dtfL2.p4() : dtfL1.p4() ;
    //lhcb_angles(bzero, dtf_eplus, dtf_eminus, dtfKaon.p4(), dtfPion.p4(), &ctl_dtf, &ctk_dtf, &phi_dtf);

    // Add recovered photons           M.Soares 20.Oct.21
    if(use_MC==true ) 
    {
      //// Run it only for correct matches!! 
      if( abs(trueB->ID) == true_ids[channel][0] && abs(trueKaon->ID) == true_ids[channel][1] && abs(truePion->ID) == true_ids[channel][2] && abs(trueL1->ID) == true_ids[channel][3] &&  abs(trueL2->ID) == true_ids[channel][4]  ){
          reco_wrongID = 0;
      }
      else  reco_wrongID = 1;

      /* --------> Commented out because the correct variables are are the MCDT  A.Biolchini & M.Soares 3.02.2023
        // // Run it only for correct matches!! 
        // if( abs(trueB->ID) == true_ids[channel][0] && abs(trueKaon->ID) == true_ids[channel][1] && abs(truePion->ID) == true_ids[channel][2] && abs(trueL1->ID) == true_ids[channel][3] &&  abs(trueL2->ID) == true_ids[channel][4]  ){
      	// +ve -> electron, -ve -> positron
      	const auto eplus_true  = trueL1->ID<0 ? trueL1->p4() : trueL2->p4() ;
      	const auto eminus_true = trueL1->ID<0 ? trueL2->p4() : trueL1->p4() ;
      	const auto bzero_true = trueKaon->ID > 0 ;

      	lhcb_angles(bzero_true, eplus_true, eminus_true, trueKaon->p4(), truePion->p4(), ctl_lhcb_true, ctk_lhcb_true, phi_lhcb_true);

      	// --- recover photons
      	//
      	const auto jPsi = trueB->p4() - trueKstar->p4() ;
      	const auto boostToJ = jPsi.BoostToCM() ;
      	vec4d const boostJPsi = boost(jPsi, boostToJ) ;
      	vec4d const boostlepPlus = boost(eplus_true, boostToJ) ;
      	vec4d const boostlepMin = boost(eminus_true, boostToJ) ;
      	vec4d const boostphotJPsi = boostJPsi - boostlepPlus - boostlepMin ;
      	//nPhotonRecov =0 ;
      	auto lepMin  = eminus_true;
      	auto lepPlus = eplus_true;
      	const auto Ephotons = boostphotJPsi.E() ;
      	if( Ephotons > 1.0 ) { // set some lower limit (1 MeV) to the correction, to prevent divisions by zero
      	  const auto Ep = boostlepPlus.E() ;
      	  const auto Em = boostlepMin.E() ;
      	  // divide the energy such that both leptons get the same energy in the (B-K*) rest frame
      	  // restruct to range (0,1) to reduce sensitivity to numerical errors
      	  const auto wp = std::min(std::max(0.5*(Ephotons + Em - Ep) /  Ephotons, 0.0), 1.0 )  ; // restrict to range (0,1)
      	  const auto wm = 1-wp ;
      	  const auto boostphotJPsi = boostJPsi - boostlepPlus - boostlepMin ;
      	  const vec4d photEp = wp * boostphotJPsi ;
      	  const vec4d photEm = wm * boostphotJPsi ;
      	  lepPlus = eplus_true + boost(photEp, -boostToJ) ;
      	  lepMin  = eminus_true + boost(photEm, -boostToJ) ;
      	  if( wp>0.005) nPhotonRecov += 1 ;
      	  if( wm>0.005) nPhotonRecov += 1 ;
      	}
      	// q2_recov = jPsi.M2()/1e06 ;
      	*L1_RECOVP = trueL1->ID<0 ? lepPlus : lepMin ;
      	*L2_RECOVP = trueL1->ID<0 ? lepMin : lepPlus ;
      	
      	lhcb_angles(bzero_true, lepPlus, lepMin, trueKaon->p4(), truePion->p4(), ctl_lhcb_recov, ctk_lhcb_recov, phi_lhcb_recov);
      
      //----------------------------------------------------------------------------------------------------------------
      */

        if (channel == "B02Kst0Jpsi2ee" || channel == "B02Kst0eeSignal")
        {
          const auto MCDT_eplus_true  = MCDT_trueL1->ID<0 ? MCDT_trueL1->p4() : MCDT_trueL2->p4() ;
          const auto MCDT_eminus_true = MCDT_trueL1->ID<0 ? MCDT_trueL2->p4() : MCDT_trueL1->p4() ;
          const auto MCDT_bzero_true = MCDT_trueKaon->ID > 0 ;

          lhcb_angles(MCDT_bzero_true, MCDT_eplus_true, MCDT_eminus_true, MCDT_trueKaon->p4(), MCDT_truePion->p4(), MCDT_ctl_lhcb_true, MCDT_ctk_lhcb_true, MCDT_phi_lhcb_true);

          // --- recover photons
          //
          const auto MCDT_jPsi = MCDT_trueB->p4() - MCDT_trueKstar->p4() ;
          const auto MCDT_boostToJ = MCDT_jPsi.BoostToCM() ;
          vec4d const MCDT_boostJPsi = boost(MCDT_jPsi, MCDT_boostToJ) ;
          vec4d const MCDT_boostlepPlus = boost(MCDT_eplus_true, MCDT_boostToJ) ;
          vec4d const MCDT_boostlepMin = boost(MCDT_eminus_true, MCDT_boostToJ) ;
          vec4d const MCDT_boostphotJPsi = MCDT_boostJPsi - MCDT_boostlepPlus - MCDT_boostlepMin ;
          MCDT_nPhotonRecov = 0 ;
          auto MCDT_lepMin  = MCDT_eminus_true;
          auto MCDT_lepPlus = MCDT_eplus_true;
          const auto MCDT_Ephotons = MCDT_boostphotJPsi.E() ;
          if( MCDT_Ephotons > 1.0 ) { // set some lower limit (1 MeV) to the correction, to prevent divisions by zero
            const auto MCDT_Ep = MCDT_boostlepPlus.E() ;
            const auto MCDT_Em = MCDT_boostlepMin.E() ;
            // divide the energy such that both leptons get the same energy in the (B-K*) rest frame
            // restruct to range (0,1) to reduce sensitivity to numerical errors
            const auto MCDT_wp = std::min(std::max(0.5*(MCDT_Ephotons + MCDT_Em - MCDT_Ep) /  MCDT_Ephotons, 0.0), 1.0 )  ; // restrict to range (0,1)
            const auto MCDT_wm = 1-MCDT_wp ;
            const auto MCDT_boostphotJPsi = MCDT_boostJPsi - MCDT_boostlepPlus - MCDT_boostlepMin ;
            const vec4d MCDT_photEp = MCDT_wp * MCDT_boostphotJPsi ;
            const vec4d MCDT_photEm = MCDT_wm * MCDT_boostphotJPsi ;
            MCDT_lepPlus = MCDT_eplus_true + boost(MCDT_photEp, -MCDT_boostToJ) ;
            MCDT_lepMin  = MCDT_eminus_true + boost(MCDT_photEm, -MCDT_boostToJ) ;
            if( MCDT_wp>0.005) MCDT_nPhotonRecov += 1 ;
            if( MCDT_wm>0.005) MCDT_nPhotonRecov += 1 ;
          }
          // q2_recov = jPsi.M2()/1e06 ;
          *MCDT_L1_RECOVP = MCDT_trueL1->ID<0 ? MCDT_lepPlus : MCDT_lepMin ;
          *MCDT_L2_RECOVP = MCDT_trueL1->ID<0 ? MCDT_lepMin : MCDT_lepPlus ;
          
          lhcb_angles(MCDT_bzero_true, MCDT_lepPlus, MCDT_lepMin, MCDT_trueKaon->p4(), MCDT_truePion->p4(), MCDT_ctl_lhcb_recov, MCDT_ctk_lhcb_recov, MCDT_phi_lhcb_recov);
        }

        // Add columns with smeared electron energy
        const double landausigma[2] = {  0.0019, 0.00250 } ;
        const double landaumu[2]    = { -0.0029, -0.00245 } ;
        const int bremcatL1 = L1->BremMultiplicity >0 ;
        const int bremcatL2 = L2->BremMultiplicity >0 ;
        const auto d1 = - (landaumu[bremcatL1] + landausigma[bremcatL1] * throwlandautruncated() ) ;
        const auto d2 = - (landaumu[bremcatL2] + landausigma[bremcatL2] * throwlandautruncated() ) ;
        const auto L1p4smeared = (1+d1) * L1->p4() ;
        const auto L2p4smeared = (1+d2) * L2->p4() ;
        const auto dileptonp4smeared = L1p4smeared + L2p4smeared ;
        const auto dileptonp4 = L1->p4() + L2->p4() ;
        const auto kstarp4 = Kstar->p4() ;

        // we implement this as delta's on the mass, such that we can exploit the improvement on the mass from the vertex constraints
        const auto DeltaB_M    = (dileptonp4smeared + kstarp4).M() - (dileptonp4 + kstarp4).M() ;
        const auto DeltaJpsi_M = dileptonp4smeared.M() - dileptonp4.M() ;
        B_M_smeared    = Double_t(B->M) + DeltaB_M ;
        Jpsi_M_smeared = Double_t(JPsi->M) + DeltaJpsi_M ;
        *L1_SMEARP = L1p4smeared;
        *L2_SMEARP = L2p4smeared;

    }
    /////////////////////////////////////////New M.Geijsen & A.Biolchini & M.Soares Dec-jan 2022-2023
    //
    
    else{

    B_M_Tr_SubstL1_e2pi  = ( Trkaon->p4(mkaon) + Trpion->p4(mpion) + TrL1->p4(mpion) + TrL2->p4(melec) ).M() ;
    B_M_Tr_SubstL1L2_ee2pipi = ( Trkaon->p4(mkaon) + Trpion->p4(mpion) + TrL1->p4(mpion) + TrL2->p4(mpion) ).M() ;
    B_M_Tr_SubstL1L2_ee2kk   = ( Trkaon->p4(mkaon) + Trpion->p4(mpion) + TrL1->p4(mkaon) + TrL2->p4(mkaon) ).M() ;

    }
    //debugging
    //std::cout<< "ctl after:"<<ctl_lhcb<<std::endl;
    //std::cout<<B_M_SubstL1_e2pi<<std::endl;
    //
    ///////////////////////////////////////END

    newtree->Fill();

  }//for Nevent

  newtree->Print();
  newtree->AutoSave();
  newfile->Close();
  return;
}

void lhcb_angles(bool bzero, const vec4d & eplus, const vec4d & eminus, const vec4d & kaon, const vec4d & pion, double& costhetal, double& costhetak, double& phi)

{
  //set up boost vectors
  vec4d const b = eplus + eminus + kaon + pion;
  vec4d const ee = eplus + eminus;
  vec4d const kpi = kaon + pion;
  vec3d const eeboost(ee.BoostToCM());
  vec3d const kpiboost(kpi.BoostToCM());
  vec3d const bboost(b.BoostToCM());

  //determine costhetal
  vec4d const eminusd = boost(eminus, eeboost);
  vec4d const eplusd = boost(eplus, eeboost);
  vec4d const bd = boost(b, eeboost);

  if (bzero)
    costhetal = Cos(Angle(eplusd.Vect(), -bd.Vect()));
  else
    costhetal = Cos(Angle(eminusd.Vect(), -bd.Vect()));
  //std::cout << "Processing data" << std::endl;


  //determine costhetak
  vec4d const kaondd = boost(kaon, kpiboost);
  vec4d const bdd = boost(b, kpiboost);
  costhetak = Cos(Angle(kaondd.Vect(), -bdd.Vect()));

  //determine phi
  vec4d const kaonddd = boost(kaon, bboost);
  vec4d const pionddd = boost(pion, bboost);
  vec4d const eminusddd = boost(eminus, bboost);
  vec4d const eplusddd = boost(eplus, bboost);
  vec3d const normalkpi = kaonddd.Vect().Cross(pionddd.Vect());
  vec3d const normalee = eplusddd.Vect().Cross(eminusddd.Vect());

  vec4d const kpiddd = boost(kpi, bboost);
  if (bzero)
    {
      phi = Angle(normalkpi,normalee);  //costheta = norm.norm, all in b frame
      if (normalee.Cross(normalkpi).Dot(kpiddd.Vect()) < 0.0)
	phi = -phi;
    }
  else
    {
      phi = Angle(normalkpi, -normalee);  //costheta = norm.norm, all in b frame
      if (normalee.Cross(normalkpi).Dot(kpiddd.Vect()) < 0.0)
        phi = -phi;
    } 
};
int main(int argc, char** argv) {
    if (argc != 6) {
        cout << "Usage: " << argv[0] << " <inputFile> <outputFile> <treeName> <flag> <someString>" << endl;
        return 1;
    }

    TString inputFile = argv[1];
    TString outputFile = argv[2];
    TString treeName = argv[3];
    bool flag = (string(argv[4]) == "true");
    string someString = argv[5];

    add_angles_e(inputFile, outputFile, treeName, flag, someString);

    return 0;
}
