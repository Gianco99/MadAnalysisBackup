#include "SampleAnalyzer/User/Analyzer/user.h"
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

vector<RecJetFormat> filter_Ak4jets(vector<RecJetFormat> jets, float ptmin, float etamax, const vector<RecLeptonFormat>& leptons, float deltaRmin, const std::string& year) {

    vector<RecJetFormat> filtered;

    struct MuSigma {
        float mu_old;
        float sigma_old;
        float mu_new;
        float sigma_new;
    };
    
    // Mu and Sigma of the Log-Normal Distributions of Jet Resolution
    std::map<std::string, std::map<std::string, MuSigma>> mu_sigma_map = {
        {"<0.5", {
            {"2018", {0.33, 0.50, 0.20, 0.29}},
            {"2017", {0.33, 0.50, 0.20, 0.28}},
            {"2016", {0.33, 0.50, 0.18, 0.26}},
            {"2016APV", {0.33, 0.50, 0.17, 0.26}}
        }},
        {"0.5-1.5", {
            {"2018", {0.31, 0.48, 0.26, 0.32}},
            {"2017", {0.31, 0.48, 0.26, 0.33}},
            {"2016", {0.31, 0.48, 0.23, 0.30}},
            {"2016APV", {0.31, 0.48, 0.22, 0.29}}
        }},
        {">1.5", {
            {"2018", {0.34, 0.54, 0.33, 0.36}},
            {"2017", {0.34, 0.54, 0.32, 0.36}},
            {"2016", {0.34, 0.54, 0.29, 0.33}},
            {"2016APV", {0.34, 0.54, 0.26, 0.32}}
        }}
    };

    for (auto & jet : jets) {
        float abs_eta = fabs(jet.eta());
        std::string eta_bin;

        if (abs_eta < 0.5) {
            eta_bin = "<0.5";
        } else if (abs_eta < 1.5) {
            eta_bin = "0.5-1.5";
        } else {
            eta_bin = ">1.5";
        }

        MuSigma ms = mu_sigma_map[eta_bin][year];

        // Small minimum to avoid instability with logs
        if (jet.pt() < 5.0) {
            continue;
        }

        std::cout << "Year: " << year << std::endl;
        std::cout << "Jet pt before smearing: " << jet.pt() << std::endl;
        //std::cout << "Jet eta before smearing: " << jet.eta() << std::endl;
        //std::cout << "Jet phi before smearing: " << jet.phi() << std::endl;
        
        
        // Smear the jet pt
        float log_pt_old = log(jet.pt());
        float pt_new = exp(ms.mu_new + (ms.sigma_new / ms.sigma_old) * (log_pt_old - ms.mu_old));

        //std::cout << "Jet pt after smearing: " << pt_new << std::endl;

        MAdouble64 pt_smeared = pt_new;  // Smeared pt value
        MAdouble64 eta = jet.eta();   // Access eta from jet object
        MAdouble64 phi = jet.phi();   // Access phi from jet object
        MAdouble64 mass = jet.m();

        MAdouble64 px = pt_smeared * std::cos(phi);
        MAdouble64 py = pt_smeared * std::sin(phi);
        MAdouble64 pz = pt_smeared * std::sinh(eta);

        MAdouble64 E = std::sqrt(pt_smeared*pt_smeared * (1 + std::sinh(eta)*std::sinh(eta)) + mass*mass);

        // Create the Lorentz vector and set its components
        MALorentzVector jet_vector;
        jet_vector.SetPxPyPzE(px, py, pz, E);

        // Set the momentum for the jet
        jet.setMomentum(jet_vector);

        std::cout << "Jet pt after fancy stuff: " << jet.pt() << std::endl;
        //std::cout << "Jet eta after fancy stuff: " << jet.eta() << std::endl;
        //std::cout << "Jet phi after fancy stuff: " << jet.phi() << std::endl;

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
                         std::vector<RecJetFormat>& Ak4jets) {

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
    
    // Sort the jet collection
    std::sort(Ak4jets.begin(), Ak4jets.end(), compareBypTJets);
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

std::pair<std::vector<fastjet::PseudoJet>, std::vector<int>> getAk15Jets(std::vector<RecTrackFormat> tracks, std::vector<RecLeptonFormat> leptons, double pt_cut, double eta_cut, double d0_cut, double dz_cut, double dr_cut) {
    std::vector<fastjet::PseudoJet> input_particles;
    for (const auto &track : tracks) 
    {
        if (track.pt() >= pt_cut && std::abs(track.eta()) <= eta_cut && std::abs(track.d0()) < d0_cut && std::abs(track.dz()) < dz_cut && track.dr(leptons.at(0)) >= dr_cut && track.dr(leptons.at(1)) >= dr_cut) 
        {
            fastjet::PseudoJet particle(track.px(), track.py(), track.pz(), track.e());
            input_particles.emplace_back(particle);
        }
    }

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.5);
    fastjet::ClusterSequence clust_seq(input_particles, jet_def);
    
    std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(0.0));

    std::vector<int> constituent_sizes;
    for (const auto &jet : inclusive_jets) 
    {
        constituent_sizes.push_back(jet.constituents().size());
    }

    return std::make_pair(inclusive_jets, constituent_sizes);
}

double calculateDeltaR(const RecJetFormat& rec_jet, const fastjet::PseudoJet& pseudo_jet) 
{
    double rec_eta = rec_jet.eta();
    double rec_phi = rec_jet.phi();

    double pseudo_eta = pseudo_jet.eta();
    double pseudo_phi = pseudo_jet.phi();

    // Map both phi values to [-pi, pi] for consistency
    if (rec_phi > M_PI) {
        rec_phi -= 2 * M_PI;
    }
    if (pseudo_phi > M_PI) {
        pseudo_phi -= 2 * M_PI;
    }

    // Calculate delta eta and delta phi
    double delta_eta = rec_eta - pseudo_eta;
    double delta_phi = std::atan2(std::sin(rec_phi - pseudo_phi), std::cos(rec_phi - pseudo_phi));

    // Debugging output
    //std::cout << "rec_eta: " << rec_eta << ", pseudo_eta: " << pseudo_eta << std::endl;
    //std::cout << "rec_phi: " << rec_phi << ", pseudo_phi: " << pseudo_phi << std::endl;
    //std::cout << "delta_eta: " << delta_eta << ", delta_phi: " << delta_phi << std::endl;

    double deltaR = std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);

    //std::cout << "deltaR: " << deltaR << std::endl;

    return deltaR;
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
  Manager()->AddCut("oneAk15Cluster");
  Manager()->AddCut("onShellZMass");
  Manager()->AddCut("minZpT");
  Manager()->AddCut("noBTag");
  Manager()->AddCut("minAk15pT");
  Manager()->AddCut("dROverlap");
  
  // ===== Histograms ===== //
  Manager()->AddHisto("ZMass", 60,0.0,300.0);
  Manager()->AddHisto("ZpT", 500,0.0,500.0);
  Manager()->AddHisto("ZETA", 40,-5.0,5.0);
  Manager()->AddHisto("ZPHI", 40,-3.14,3.14);

  Manager()->AddHisto("LeadLepPT", 300,0.0,300.0);
  Manager()->AddHisto("SubleadLepPT", 300,0.0,300.0);
  Manager()->AddHisto("LeadLepETA", 20,-2.5,2.5);
  Manager()->AddHisto("SubleadLepETA", 20,-2.5,2.5);
  Manager()->AddHisto("LeadLepPHI", 20,-3.14,3.14);
  Manager()->AddHisto("SubleadLepPHI", 20,-3.14,3.14);

  Manager()->AddHisto("LeadMuonPT", 300,0.0,300.0);
  Manager()->AddHisto("SubleadMuonPT", 300,0.0,300.0);
  Manager()->AddHisto("LeadMuonETA", 20,-2.5,2.5);
  Manager()->AddHisto("SubleadMuonETA", 20,-2.5,2.5);
  Manager()->AddHisto("LeadMuonPHI", 20,-3.14,3.14);
  Manager()->AddHisto("SubleadMuonPHI", 20,-3.14,3.14);

  Manager()->AddHisto("LeadElecPT", 300,0.0,300.0);
  Manager()->AddHisto("SubleadElecPT", 300,0.0,300.0);
  Manager()->AddHisto("LeadElecETA", 20,-2.5,2.5);
  Manager()->AddHisto("SubleadElecETA", 20,-2.5,2.5);
  Manager()->AddHisto("LeadElecPHI", 20,-3.14,3.14);
  Manager()->AddHisto("SubleadElecPHI", 20,-3.14,3.14);

  Manager()->AddHisto("NJets", 10,0.0,10.0);
  Manager()->AddHisto("LeadAk4PT", 500,0.0,500.0);
  Manager()->AddHisto("LeadAk4ETA", 40,-5.0,5.0);
  Manager()->AddHisto("LeadAk4PHI", 40,-3.14,3.14);
  Manager()->AddHisto("LeadAk4NTRACKS", 100,0.0,100.0);

  Manager()->AddHisto("LeadAk15PT", 500,0.0,500.0);
  Manager()->AddHisto("LeadAk15NTRACKS", 100,0.0,100.0);
  Manager()->AddHisto("LeadAk15ETA", 40,-5.0,5.0);
  Manager()->AddHisto("LeadAk15PHI", 40,-3.14,3.14);
  Manager()->AddHisto("LeadAk15Mass", 400,0,400);

  Manager()->AddHisto("ABCD_A", 100,21.0,121.0);
  Manager()->AddHisto("ABCD_B1", 100,21.0,121.0);
  Manager()->AddHisto("ABCD_B2", 100,21.0,121.0);
  Manager()->AddHisto("ABCD_C1", 1,14.0,21.0);
  Manager()->AddHisto("ABCD_C2", 1,0.0,14.0);
  Manager()->AddHisto("ABCD_D1", 1,14.0,21.0);
  Manager()->AddHisto("ABCD_D2", 1,0.0,14.0);
  Manager()->AddHisto("ABCD_E1", 1,14.0,21.0);
  Manager()->AddHisto("ABCD_E2", 1,0.0,14.0);

  // No problem during initialization
  return true;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool user::Execute(SampleFormat& sample, const EventFormat& event)
{
  // Object Cut definitions
  float const CUT_ELECTRON_PT_MIN  = 10;
  float const CUT_ELECTRON_ETA_MAX = 2.4;

  float const CUT_MUON_PT_MIN      = 10;
  float const CUT_MUON_ETA_MAX     = 2.5;
  float const CUT_MUON_D0          = 0.02;
  float const CUT_MUON_DZ          = 0.1;

  float const CUT_AK4JET_PT_MIN  = 30.0;
  float const CUT_AK4JET_ETA_MAX = 2.5;
  float const CUT_AK4JET_LEPTON_DELTAR = 0.4;

  // Event weight
  double weight = 1.0;
  if (!Configuration().IsNoEventWeight() && event.mc() != 0) {
      weight = event.mc()->weight();
  }

  // Random number generator for year assignment
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0, 1);

  // Year probabilities
  double p2018 = 59.9 / 137.8;
  double p2017 = 41.6 / 137.8;
  double p2016 = 16.4 / 137.8;
  double p2016APV = 19.9 / 137.8;

  // Assign year based on random number
  double rand = dis(gen);
  std::string year;
  if (rand < p2018) {
      year = "2018";
      weight *= p2018;  // Multiply by 2018 probability
  } else if (rand < p2018 + p2017) {
      year = "2017";
      weight *= p2017;  // Multiply by 2017 probability
  } else if (rand < p2018 + p2017 + p2016) {
      year = "2016";
      weight *= p2016;  // Multiply by 2016 probability
  } else {
      year = "2016APV";
      weight *= p2016APV;  // Multiply by 2016APV probability
  }

  Manager()->InitializeForNewEvent(weight);
  if (event.rec() == 0) { return true; }

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

  // Applying object selections to jet collection
  vector<RecJetFormat> Ak4jets         = filter_Ak4jets(event.rec()->jets(), CUT_AK4JET_PT_MIN, CUT_AK4JET_ETA_MAX, leptons, CUT_AK4JET_LEPTON_DELTAR, year);

  // Sort all collections before event selection
  sortCollectionsBypT(electrons, posElectrons, negElectrons, muons, posMuons, negMuons, leptons, posLeptons, negLeptons, Ak4jets);

  //Event level Selections

  // Event Cut definitions
  float const LEAD_LEPTON_PT  = 25;

  float const ZMASS_LOW = 60.0;
  float const ZMASS_HIGH = 120.0;

  float const ZPT_LOW = 25.0;

  float const TRACK_PT_MIN = 1.0;
  float const TRACK_ETA_MAX = 2.5;
  float const TRACK_D0_MAX = 0.05;
  float const TRACK_DZ_MAX = 0.05;
  float const TRACK_DR_MAX = 0.4;

  float const AK15JET_PT_MIN  = 60.0;

  float const DR_MAX  = 1.5;

  // Two OSSF Lepton + Lead Lepton pT Selection
  bool twoOSleptons = (posLeptons.size() == 1 && negLeptons.size() == 1);                  // Require exactly two opposite-sign leptons (either muons or electrons)
  bool twoSFleptons = (muons.size() == 2 || electrons.size() == 2);                        // Require exactly two muons or exactly two electrons
  bool LeadpTleptons = leptons.size() > 0 && leptons.at(0).pt() >= LEAD_LEPTON_PT;     // Require the leading lepton pT to be >= 25 GeV
  bool twoOSSFLeptons = twoOSleptons && twoSFleptons && LeadpTleptons ;                    // Concatenating cuts
  if (not Manager()->ApplyCut(twoOSSFLeptons, "twoOSSFLeptons")) return true;

  // Do Ak15 clustering 
  auto Ak15result = getAk15Jets(event.rec()->EFlowTracks(), leptons, TRACK_PT_MIN, TRACK_ETA_MAX, TRACK_D0_MAX, TRACK_DZ_MAX, TRACK_DR_MAX);
  std::vector<fastjet::PseudoJet> Ak15jets = Ak15result.first;
  std::vector<int> Ak15jetConstituentSizes = Ak15result.second;

  // One Ak15 cluster
  bool oneAk15Cluster = (Ak15jets.size() > 0);
  if (not Manager()->ApplyCut(oneAk15Cluster, "oneAk15Cluster")) return true;

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

  bool dROverlap = false;
  if (!Ak4jets.empty()) {
      dROverlap = (calculateDeltaR(Ak4jets.at(0), Ak15jets.at(0)) < DR_MAX);
  }
  //std::cout << "dROverlap value: " << dROverlap << std::endl;
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

  Manager()->FillHisto("NJets", Ak4jets.size());
  Manager()->FillHisto("LeadAk4PT", Ak4jets.at(0).pt());
  Manager()->FillHisto("LeadAk4ETA", Ak4jets.at(0).eta());
  Manager()->FillHisto("LeadAk4PHI", Ak4jets.at(0).phi());
  Manager()->FillHisto("LeadAk4NTRACKS", Ak4jets.at(0).ntracks());

  Manager()->FillHisto("LeadAk15PT", Ak15jets.at(0).pt());
  Manager()->FillHisto("LeadAk15NTRACKS", Ak15jetConstituentSizes.at(0));
  Manager()->FillHisto("LeadAk15ETA", Ak15jets.at(0).eta());
  double phi_mapped = Ak15jets.at(0).phi();
  if (phi_mapped > M_PI) {
      phi_mapped -= 2 * M_PI;
  }
  Manager()->FillHisto("LeadAk15PHI", phi_mapped);
  Manager()->FillHisto("LeadAk15Mass", Ak15jets.at(0).m());

  if( (Ak15jetConstituentSizes.at(0) >= 21.0) && (Ak4jets.at(0).pt() <= 135.0)){ Manager()->FillHisto("ABCD_A", Ak15jetConstituentSizes.at(0));}
  if( (Ak15jetConstituentSizes.at(0) >= 14.0) && (Ak15jetConstituentSizes.at(0) < 21.0) && (Ak4jets.at(0).pt() <= 135.0)){ Manager()->FillHisto("ABCD_C1", Ak15jetConstituentSizes.at(0));}
  if( (Ak15jetConstituentSizes.at(0) < 14.0) && (Ak4jets.at(0).pt() <= 135.0)){ Manager()->FillHisto("ABCD_C2", Ak15jetConstituentSizes.at(0));}

  if( (Ak15jetConstituentSizes.at(0) >= 21.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0) ){ Manager()->FillHisto("ABCD_B1", Ak15jetConstituentSizes.at(0));}
  if( (Ak15jetConstituentSizes.at(0) >= 14.0) && (Ak15jetConstituentSizes.at(0) < 21.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0) ){ Manager()->FillHisto("ABCD_D1", Ak15jetConstituentSizes.at(0));}
  if( (Ak15jetConstituentSizes.at(0) < 14.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0) ){ Manager()->FillHisto("ABCD_D2", Ak15jetConstituentSizes.at(0));}

  if( (Ak15jetConstituentSizes.at(0) >= 21.0) && (Ak4jets.at(0).pt() > 220.0) ){ Manager()->FillHisto("ABCD_B2", Ak15jetConstituentSizes.at(0));}
  if( (Ak15jetConstituentSizes.at(0) >= 14.0) && (Ak15jetConstituentSizes.at(0) < 21.0) && (Ak4jets.at(0).pt() > 220.0) ){ Manager()->FillHisto("ABCD_E1", Ak15jetConstituentSizes.at(0));}
  if( (Ak15jetConstituentSizes.at(0) < 14.0)  && (Ak4jets.at(0).pt() > 220.0) ){ Manager()->FillHisto("ABCD_E2", Ak15jetConstituentSizes.at(0));}

  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void user::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
}


