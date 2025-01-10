#include "SampleAnalyzer/User/Analyzer/cms_exo_23_003.h"
using namespace MA5;
using namespace std;

vector<RecLeptonFormat> filter_muons(vector<RecLeptonFormat> objects, float ptmin, float etamax, float d0, float dz, std::string charge = "") {
  // Helper function to select muons passing pt, eta, and impact parameter selections
  vector<RecLeptonFormat> filtered;

  for (auto & obj : objects) {

    // Charge filter
    if (charge == "+") {
      if (obj.charge() <= 0) {
        continue; // Only keep positive charge if "+" is specified
      }
    } else if (charge == "-") {
      if (obj.charge() >= 0) {
        continue; // Only keep negative charge if "-" is specified
      }
    }

    // Object Selections
    if (obj.pt() < ptmin) { // pT
      continue;
    }
    if (fabs(obj.eta()) > etamax) { // eta
      continue;
    }
    if (fabs(obj.d0()) > d0) { //d0
      continue;
    }
    if (fabs(obj.dz()) > dz) { //dz
      continue;
    }
    filtered.push_back(obj);
  }

  return filtered;
}

vector<RecLeptonFormat> filter_electrons(vector<RecLeptonFormat> objects, float ptmin, float etamax, std::string charge = "") {
  // Helper function to select electrons passing pt, eta, and impact parameter selections
  vector<RecLeptonFormat> filtered;

  for (auto & obj : objects) {

    // Charge filter
    if (charge == "+") {
      if (obj.charge() <= 0) {
        continue; // Only keep positive charge if "+" is specified
      }
    } else if (charge == "-") {
      if (obj.charge() >= 0) {
        continue; // Only keep negative charge if "-" is specified
      }
    }

    // Object Selections
    double abseta = fabs(obj.eta());

    if (obj.pt() < ptmin) {
      continue;
    }
    if (abseta > etamax) {
      continue;
    }

    // Electron-specific eta and impact parameter cuts
    // Reject electrons with ABSETA > 1.444 and ABSETA < 1.566
    if (abseta > 1.444 && abseta < 1.566) {
      continue;
    }

    // Apply impact parameter cuts based on eta region
    if(fabs(obj.d0()) > (0.05 + 0.05*(abseta > 1.479))){
      continue;
    }

    if(fabs(obj.dz()) > (0.1 + 0.1*(abseta > 1.479))){
      continue;
    }

    filtered.push_back(obj);
  }

  return filtered;
}

vector<RecJetFormat> filter_Ak4jets(vector<RecJetFormat> jets, float ptmin, float etamax, const vector<RecLeptonFormat>& leptons, float deltaRmin) {
  // Helper function to select jets passing pt, eta, and ΔR selections
  vector<RecJetFormat> filtered;
  
  for (auto & jet : jets) {

    // Object Selections
    if (jet.pt() < ptmin) {
      continue;
    }
    if (fabs(jet.eta()) > etamax) {
      continue;
    }
    
    // ΔR selection (jet-lepton separation)
    bool passesDeltaRCut = true;
    
    // Loop over each lepton to check if any are too close to the jet
    for (const auto & lepton : leptons) {
      // Calculate ΔR between the jet and the lepton
      if (jet.dr(lepton) <= deltaRmin) {
        passesDeltaRCut = false;
        break; // If one lepton is too close, reject this jet
      }
    }
    
    // If the jet passes the ΔR cut, add it to the filtered list
    if (passesDeltaRCut) {
      filtered.push_back(jet);
    }
  }

  return filtered;
}

vector<RecJetFormat> filter_Ak15jets(vector<RecJetFormat> objects, float ptmin, float etamax) {
  // Helper function to select jets passing pt and eta selections
  vector<RecJetFormat> filtered;
  
  for (auto & obj : objects) {
    
    // Object Selections
    if (obj.pt() < ptmin) { // 
      continue;
    }
    if (fabs(obj.eta()) > etamax) {
      continue;
    }
    
    filtered.push_back(obj);
  }

  return filtered;
}

bool compareBypT(const RecLeptonFormat &a, const RecLeptonFormat &b) {
    return a.pt() > b.pt();  // Sort by pT in descending order
}

bool compareBypTJets(const RecJetFormat &a, const RecJetFormat &b) {
    return a.pt() > b.pt();  // Sort by pT in descending order
}

void sortCollectionsBypT(std::vector<RecLeptonFormat>& electrons, 
                         std::vector<RecLeptonFormat>& posElectrons, 
                         std::vector<RecLeptonFormat>& negElectrons, 
                         std::vector<RecLeptonFormat>& muons, 
                         std::vector<RecLeptonFormat>& posMuons, 
                         std::vector<RecLeptonFormat>& negMuons, 
                         std::vector<RecLeptonFormat>& leptons, 
                         std::vector<RecLeptonFormat>& posLeptons, 
                         std::vector<RecLeptonFormat>& negLeptons, 
                         std::vector<RecJetFormat>& Ak4jets, 
                         std::vector<RecJetFormat>& Ak15jets) {

    // Sort the electron and muon collections
    std::sort(electrons.begin(), electrons.end(), compareBypT);
    std::sort(posElectrons.begin(), posElectrons.end(), compareBypT);
    std::sort(negElectrons.begin(), negElectrons.end(), compareBypT);
    
    std::sort(muons.begin(), muons.end(), compareBypT);
    std::sort(posMuons.begin(), posMuons.end(), compareBypT);
    std::sort(negMuons.begin(), negMuons.end(), compareBypT);
    
    // Sort the combined lepton collections
    std::sort(leptons.begin(), leptons.end(), compareBypT);
    std::sort(posLeptons.begin(), posLeptons.end(), compareBypT);
    std::sort(negLeptons.begin(), negLeptons.end(), compareBypT);
    
    // Sort the jet collections
    std::sort(Ak4jets.begin(), Ak4jets.end(), compareBypTJets);
    std::sort(Ak15jets.begin(), Ak15jets.end(), compareBypTJets);
}

bool BTagVeto(vector<RecJetFormat> jets){
  // Returns true if a b jet is found
  for(auto jet : jets) {
    if (jet.btag()) {
      return false;
    }
  }
  return true;
}

// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
MAbool user::Initialize(const MA5::Configuration& cfg,
                      const std::map<std::string,std::string>& parameters)
{
  // Initializing PhysicsService for RECO
  PHYSICS->recConfig().Reset();

  PHYSICS->recConfig().UseDeltaRIsolation(0.5);

  // ===== Signal region ===== //
  Manager()->AddRegionSelection("SR");

  // ===== Selections ===== //
  Manager()->AddCut("twoOSSFLeptons");
  Manager()->AddCut("oneOfEachJet");
  Manager()->AddCut("onShellZMass");
  Manager()->AddCut("minZpT");
  Manager()->AddCut("noBTag");
  Manager()->AddCut("minAk15pT");
  Manager()->AddCut("dROverlap");
  
  // ===== Histograms ===== //
  Manager()->AddHisto("ZMass", 60,60.0,120.0);
  Manager()->AddHisto("ZpT", 500,0.0,500.0);
  Manager()->AddHisto("ZETA", 60,-3.0,3.0);
  Manager()->AddHisto("ZPHI", 60,-3.14,3.14);

  Manager()->AddHisto("LeadLepPT", 300,0.0,300.0);
  Manager()->AddHisto("SubleadLepPT", 300,0.0,300.0);
  Manager()->AddHisto("LeadLepETA", 60,-3.0,3.0);
  Manager()->AddHisto("SubleadLepETA", 60,-3.0,3.0);
  Manager()->AddHisto("LeadLepPHI", 60,-3.14,3.14);
  Manager()->AddHisto("SubleadLepPHI", 60,-3.14,3.14);

  Manager()->AddHisto("LeadMuonPT", 300,0.0,300.0);
  Manager()->AddHisto("SubleadMuonPT", 300,0.0,300.0);
  Manager()->AddHisto("LeadMuonETA", 60,-3.0,3.0);
  Manager()->AddHisto("SubleadMuonETA", 60,-3.0,3.0);
  Manager()->AddHisto("LeadMuonPHI", 60,-3.14,3.14);
  Manager()->AddHisto("SubleadMuonPHI", 60,-3.14,3.14);

  Manager()->AddHisto("LeadElecPT", 300,0.0,300.0);
  Manager()->AddHisto("SubleadElecPT", 300,0.0,300.0);
  Manager()->AddHisto("LeadElecETA", 60,-3.0,3.0);
  Manager()->AddHisto("SubleadElecETA", 60,-3.0,3.0);
  Manager()->AddHisto("LeadElecPHI", 60,-3.14,3.14);
  Manager()->AddHisto("SubleadElecPHI", 60,-3.14,3.14);

  Manager()->AddHisto("LeadAk4PT", 500,0.0,500.0);
  Manager()->AddHisto("LeadAk4ETA", 60,-3.0,3.0);
  Manager()->AddHisto("LeadAk4PHI", 60,-3.14,3.14);
  Manager()->AddHisto("LeadAk4NTRACKS", 100,0.0,100.0);

  Manager()->AddHisto("LeadAk15PT", 500,0.0,500.0);
  Manager()->AddHisto("LeadAk15NTRACKS", 100,0.0,100.0);
  Manager()->AddHisto("LeadAk15ETA", 60,-3.0,3.0);
  Manager()->AddHisto("LeadAk15PHI", 60,-3.14,3.14);
  Manager()->AddHisto("LeadAk15Mass", 400,0,400);



  // No problem during initialization
  return true;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool cms_exo_23_003::Execute(SampleFormat& sample, const EventFormat& event)
{
  // Object Cut definitions
  float const CUT_ELECTRON_PT_MIN  = 10;
  float const CUT_ELECTRON_ETA_MAX = 2.4;

  float const CUT_MUON_PT_MIN      = 10;
  float const CUT_MUON_ETA_MAX     = 2.5;
  float const CUT_MUON_D0          = 0.02;
  float const CUT_MUON_DZ          = 0.1;

  float const CUT_AK4JET_PT_MIN  = 20.0;
  float const CUT_AK4JET_ETA_MAX = 2.5;
  float const CUT_AK4JET_LEPTON_DELTAR = 0.4;

  // Event weight
  double weight=1.;
  if (!Configuration().IsNoEventWeight() && event.mc()!=0) {
    weight=event.mc()->weight();
  }
  Manager()->InitializeForNewEvent(weight);
  if (event.rec()==0) {return true;}

  // Applying object selections to lepton collections
  vector<RecLeptonFormat> electrons    = filter_electrons(event.rec()->electrons(), CUT_ELECTRON_PT_MIN, CUT_ELECTRON_ETA_MAX);
  vector<RecLeptonFormat> posElectrons = filter_electrons(event.rec()->electrons(), CUT_ELECTRON_PT_MIN, CUT_ELECTRON_ETA_MAX, "+");
  vector<RecLeptonFormat> negElectrons = filter_electrons(event.rec()->electrons(), CUT_ELECTRON_PT_MIN, CUT_ELECTRON_ETA_MAX, "-");

  vector<RecLeptonFormat> muons        = filter_muons(event.rec()->muons(), CUT_MUON_PT_MIN, CUT_MUON_ETA_MAX, CUT_MUON_D0, CUT_MUON_DZ);
  vector<RecLeptonFormat> posMuons     = filter_muons(event.rec()->muons(), CUT_MUON_PT_MIN, CUT_MUON_ETA_MAX, CUT_MUON_D0, CUT_MUON_DZ, "+");
  vector<RecLeptonFormat> negMuons     = filter_muons(event.rec()->muons(), CUT_MUON_PT_MIN, CUT_MUON_ETA_MAX, CUT_MUON_D0, CUT_MUON_DZ, "-");

  vector<RecLeptonFormat> leptons = electrons; // Start with electrons
  leptons.insert(leptons.end(), muons.begin(), muons.end()); // Add muons
  
  vector<RecLeptonFormat> posLeptons = posElectrons; // Start with positive electrons
  posLeptons.insert(posLeptons.end(), posMuons.begin(), posMuons.end()); // Add positive muons
  
  vector<RecLeptonFormat> negLeptons = negElectrons; // Start with negative electrons
  negLeptons.insert(negLeptons.end(), negMuons.begin(), negMuons.end()); // Add negative muons

  // Applying object selections to jet collections 
  vector<RecJetFormat> Ak4jets         = filter_Ak4jets(event.rec()->jets(), CUT_AK4JET_PT_MIN, CUT_AK4JET_ETA_MAX, leptons, CUT_AK4JET_LEPTON_DELTAR);
  vector<RecJetFormat> Ak15jets        = filter_Ak15jets(event.rec()->fatjets(), 0.0, 999);

  // Sort all collections before event selection
  sortCollectionsBypT(electrons, posElectrons, negElectrons, muons, posMuons, negMuons, leptons, posLeptons, negLeptons, Ak4jets, Ak15jets);

  //Event level Selections

  // Event Cut definitions
  float const LEAD_LEPTON_PT  = 25;

  float const ZMASS_LOW = 60.0;
  float const ZMASS_HIGH = 120.0;

  float const ZPT_LOW = 25.0;

  float const AK15JET_PT_MIN  = 60.0;

  float const DR_MAX  = 1.5;

  // Two OSSF Lepton + Lead Lepton pT Selection
  bool twoOSleptons = (posLeptons.size() == 1 && negLeptons.size() == 1);                  // Require exactly two opposite-sign leptons (either muons or electrons)
  bool twoSFleptons = (muons.size() == 2 || electrons.size() == 2);                        // Require exactly two muons or exactly two electrons
  bool LeadpTleptons = leptons.size() > 0 && leptons.at(0).pt() >= LEAD_LEPTON_PT;     // Require the leading lepton pT to be >= 25 GeV
  bool twoOSSFLeptons = twoOSleptons && twoSFleptons && LeadpTleptons ;                    // Concatenating cuts
  if (not Manager()->ApplyCut(twoOSSFLeptons, "twoOSSFLeptons")) return true;

  // One Ak15 cluster and One Ak4 jet Selection 
  bool oneAk15Cluster = (Ak15jets.size() > 0);
  bool oneAk4Cluster = (Ak4jets.size() > 0);
  bool oneOfEachJet = (oneAk15Cluster && oneAk4Cluster);
  if (not Manager()->ApplyCut(oneOfEachJet, "oneOfEachJet")) return true;

  // Reconstruct the Z

  // At this point there are exactly one positive and one negative lepton
  RecLeptonFormat posLepton = posLeptons.at(0);
  RecLeptonFormat negLepton = negLeptons.at(0);

  ParticleBaseFormat recoZ;
  recoZ += posLepton.momentum();
  recoZ += negLepton.momentum();

  // On-Shell ZMass Selection
  bool onShellZMass = (recoZ.m() >= ZMASS_LOW && recoZ.m() <= ZMASS_HIGH);
  if (not Manager()->ApplyCut(onShellZMass, "onShellZMass")) return true;

  // ZpT Selection
  bool minZpT = (recoZ.pt() >= ZPT_LOW);
  if (not Manager()->ApplyCut(minZpT, "minZpT")) return true;

  // BTag Veto
  bool noBTag = BTagVeto(Ak4jets);
  if (not Manager()->ApplyCut(noBTag, "noBTag")) return true;

  // Ak15 pT
  bool minAk15pT = (Ak15jets.at(0).pt() > AK15JET_PT_MIN);
  if (not Manager()->ApplyCut(minAk15pT, "minAk15pT")) return true;

  // dR overlap between lead Ak4 and lead Ak15
  bool dROverlap = (Ak4jets.at(0).dr(Ak15jets.at(0)) < DR_MAX);
  if (not Manager()->ApplyCut(dROverlap, "dROverlap")) return true;

  // Fill histograms
  Manager()->FillHisto("ZMass", recoZ.m());
  Manager()->FillHisto("ZpT",  recoZ.pt());
  Manager()->FillHisto("ZETA", recoZ.eta());
  Manager()->FillHisto("ZPHI",  recoZ.phi());

  Manager()->FillHisto("LeadLepPT", leptons.at(0).pt());
  Manager()->FillHisto("SubleadLepPT", leptons.at(1).pt());
  Manager()->FillHisto("LeadLepETA", leptons.at(0).eta());
  Manager()->FillHisto("SubleadLepETA", leptons.at(1).eta());
  Manager()->FillHisto("LeadLepPHI", leptons.at(0).phi());
  Manager()->FillHisto("SubleadLepPHI", leptons.at(1).phi());

  if(muons.size() > 0){
    Manager()->FillHisto("LeadMuonPT", muons.at(0).pt());
    Manager()->FillHisto("SubleadMuonPT", muons.at(1).pt());
    Manager()->FillHisto("LeadMuonETA", muons.at(0).eta());
    Manager()->FillHisto("SubleadMuonETA", muons.at(1).eta());
    Manager()->FillHisto("LeadMuonPHI", muons.at(0).phi());
    Manager()->FillHisto("SubleadMuonPHI", muons.at(1).phi());
  }
  
  else{
    Manager()->FillHisto("LeadElecPT", electrons.at(0).pt());
    Manager()->FillHisto("SubleadElecPT", electrons.at(1).pt());
    Manager()->FillHisto("LeadElecETA", electrons.at(0).eta());
    Manager()->FillHisto("SubleadElecETA", electrons.at(1).eta());
    Manager()->FillHisto("LeadElecPHI", electrons.at(0).phi());
    Manager()->FillHisto("SubleadElecPHI", electrons.at(1).phi());
  }

  Manager()->FillHisto("LeadAk4PT", Ak4jets.at(0).pt());
  Manager()->FillHisto("LeadAk4ETA", Ak4jets.at(0).eta());
  Manager()->FillHisto("LeadAk4PHI", Ak4jets.at(0).phi());
  Manager()->FillHisto("LeadAk4NTRACKS", Ak4jets.at(0).ntracks());

  Manager()->FillHisto("LeadAk15PT", Ak15jets.at(0).pt());
  Manager()->FillHisto("LeadAk15NTRACKS", Ak15jets.at(0).ntracks());
  Manager()->FillHisto("LeadAk15ETA", Ak15jets.at(0).eta());
  Manager()->FillHisto("LeadAk15PHI", Ak15jets.at(0).phi());
  Manager()->FillHisto("LeadAk15Mass", Ak15jets.at(0).m());

  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void cms_exo_23_003::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
}


