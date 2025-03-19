#include "SampleAnalyzer/User/Analyzer/user.h"
using namespace MA5;
using namespace std;
#include <random>
#include <Eigen/Dense>

double sphericity(const std::vector<fastjet::PseudoJet>& particles, double r) {
    if (particles.empty()) {
        std::cerr << "No particles in the leading jet!" << std::endl;
        return 0.0;
    }

    // Initialize sums for the sphericity matrix
    double S_xx = 0.0, S_xy = 0.0, S_xz = 0.0;
    double S_yy = 0.0, S_yz = 0.0, S_zz = 0.0;
    double norm = 0.0;

    // Calculate momentum components and normalization factor
    for (const auto& particle : particles) {
        double px = particle.px();
        double py = particle.py();
        double pz = particle.pz();
        double p = std::sqrt(px*px + py*py + pz*pz);  // Magnitude of momentum

        if (p == 0) continue;

        double weight = std::pow(p, r - 2.0);  // Weight for each component

        // Sphericity tensor components
        S_xx += px * px * weight;
        S_xy += px * py * weight;
        S_xz += px * pz * weight;
        S_yy += py * py * weight;
        S_yz += py * pz * weight;
        S_zz += pz * pz * weight;

        // Normalization
        norm += std::pow(p, r);
    }

    // Normalize sphericity matrix components
    if (norm > 0) {
        S_xx /= norm;
        S_xy /= norm;
        S_xz /= norm;
        S_yy /= norm;
        S_yz /= norm;
        S_zz /= norm;
    }

    // Construct the sphericity matrix
    Eigen::Matrix3d S;
    S << S_xx, S_xy, S_xz,
         S_xy, S_yy, S_yz,
         S_xz, S_yz, S_zz;

    // Compute eigenvalues of the sphericity matrix
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(S);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigenvalue computation failed!" << std::endl;
        return 0.0;
    }

    // Get sorted eigenvalues
    Eigen::Vector3d evals = solver.eigenvalues();
    std::sort(evals.data(), evals.data() + evals.size());

    // Extract the smallest two eigenvalues: eval1 and eval2
    double eval1 = evals[0];
    double eval2 = evals[1];

    // Calculate the sphericity value: 1.5 * (eval1 + eval2)
    double sphericity_value = 1.5 * (eval1 + eval2);

    return sphericity_value;
}

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

std::pair<std::vector<fastjet::PseudoJet>, std::vector<std::vector<fastjet::PseudoJet>>> getAk15Jets(std::vector<RecTrackFormat> tracks, std::vector<RecLeptonFormat> leptons, double pt_cut, double eta_cut, double d0_cut, double dz_cut, double dr_cut) 
{
    std::vector<fastjet::PseudoJet> input_particles;
    std::random_device rd; // Initialize random number generator
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0.0, 1.0); // Uniform distribution between 0 and 1
    
    for (const auto &track : tracks) 
    {
        // Generate a random number and check if it's less than 0.90
       // if (dis(gen) > 0.90) {
        //    continue; // Skip the track with 90% chance
        //}

        // Check if the track passes all the selection cuts
        if (track.pt() >= pt_cut && std::abs(track.eta()) <= eta_cut && 
            std::abs(track.d0()) < d0_cut && std::abs(track.dz()) < dz_cut && 
            track.dr(leptons.at(0)) >= dr_cut && track.dr(leptons.at(1)) >= dr_cut) 
        {
            fastjet::PseudoJet particle(track.px(), track.py(), track.pz(), track.e());
            input_particles.emplace_back(particle);
        }
    }

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 1.5);
    fastjet::ClusterSequence clust_seq(input_particles, jet_def);
    
    std::vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(0.0));

    std::vector<std::vector<fastjet::PseudoJet>> cluster_constituents;

    for (const auto &jet : inclusive_jets) 
    {
        cluster_constituents.push_back(jet.constituents());

    }

    return std::make_pair(inclusive_jets, cluster_constituents);
}

double calculateDeltaR(const RecJetFormat& rec_jet, const fastjet::PseudoJet& pseudo_jet) 
{
    double rec_eta = rec_jet.eta();
    double rec_phi = rec_jet.phi();

    double pseudo_eta = pseudo_jet.eta();
    double pseudo_phi = pseudo_jet.phi();

    // Map pseudo_phi from [0, 2pi] to [-pi, pi] for consistency
    if (pseudo_phi > M_PI) {
        pseudo_phi -= 2 * M_PI;
    }

    // Calculate delta eta and delta phi
    double delta_eta = rec_eta - pseudo_eta;
    double delta_phi = std::atan2(std::sin(rec_phi - pseudo_phi), std::cos(rec_phi - pseudo_phi));

    return std::sqrt(delta_eta * delta_eta + delta_phi * delta_phi);
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
  Manager()->AddRegionSelection("CRDY");

  // ===== Selections ===== //
  Manager()->AddCut("twoOSSFLeptons");
  Manager()->AddCut("oneAk15Cluster");
  Manager()->AddCut("offShellZMass");
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
  Manager()->AddHisto("LeadAk4ETA", 100,-5.0,5.0);
  Manager()->AddHisto("LeadAk4PHI", 40,-3.14,3.14);
  Manager()->AddHisto("LeadAk4NTRACKS", 100,0.0,100.0);

  Manager()->AddHisto("LeadAk15PT", 500,0.0,500.0);
  Manager()->AddHisto("LeadAk15NTRACKS", 200,0.0,200.0);
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

  Manager()->AddHisto("MET_PT", 500,0.0,500.0);
  Manager()->AddHisto("MET_ETA", 40,-5.0,5.0);
  Manager()->AddHisto("MET_PHI", 40,-3.14,3.14);

  Manager()->AddHisto("Sphericity", 50, 0.0, 1.0);

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

  // Applying object selections to jet collection
  vector<RecJetFormat> Ak4jets         = filter_Ak4jets(event.rec()->jets(), CUT_AK4JET_PT_MIN, CUT_AK4JET_ETA_MAX, leptons, CUT_AK4JET_LEPTON_DELTAR);

  // Sort all collections before event selection
  sortCollectionsBypT(electrons, posElectrons, negElectrons, muons, posMuons, negMuons, leptons, posLeptons, negLeptons, Ak4jets);

  //Event level Selections

  // Event Cut definitions
  float const LEAD_LEPTON_PT  = 25;

  float const ZMASS_LOW = 120.0;

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
  auto Ak15result = getAk15Jets(event.rec()->tracks(), leptons, TRACK_PT_MIN, TRACK_ETA_MAX, TRACK_D0_MAX, TRACK_DZ_MAX, TRACK_DR_MAX);

  std::vector<fastjet::PseudoJet> Ak15jets = Ak15result.first;
  std::vector<std::vector<fastjet::PseudoJet>> Ak15jetConstituents = Ak15result.second;

  // One Ak15 cluster
  bool oneAk15Cluster = (Ak15jets.size() > 0);
  if (not Manager()->ApplyCut(oneAk15Cluster, "oneAk15Cluster")) return true;

  double sphericity_value = sphericity(Ak15jetConstituents.at(0), 2.0);

  // Reconstruct the Z

  // At this point there are exactly one positive and one negative lepton
  RecLeptonFormat posLepton = posLeptons.at(0);
  RecLeptonFormat negLepton = negLeptons.at(0);

  ParticleBaseFormat recoZ;
  recoZ += posLepton.momentum();
  recoZ += negLepton.momentum();

  // Off-Shell ZMass Selection
  bool offShellZMass = (recoZ.m() >= ZMASS_LOW);
  if (not Manager()->ApplyCut(offShellZMass, "offShellZMass")) return true;

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
  Manager()->FillHisto("LeadAk15NTRACKS", Ak15jetConstituents.at(0).size());
  Manager()->FillHisto("LeadAk15ETA", Ak15jets.at(0).eta());
  double phi_mapped = Ak15jets.at(0).phi();
  if (phi_mapped > M_PI) {
      phi_mapped -= 2 * M_PI;
  }
  Manager()->FillHisto("LeadAk15PHI", phi_mapped);
  Manager()->FillHisto("LeadAk15Mass", Ak15jets.at(0).m());

  if( (Ak15jetConstituents.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() <= 135.0)){ Manager()->FillHisto("ABCD_A", Ak15jetConstituents.at(0).size());}
  if( (Ak15jetConstituents.at(0).size() >= 14.0) && (Ak15jetConstituents.at(0).size() < 21.0) && (Ak4jets.at(0).pt() <= 135.0)){ Manager()->FillHisto("ABCD_C1", Ak15jetConstituents.at(0).size());}
  if( (Ak15jetConstituents.at(0).size() < 14.0) && (Ak4jets.at(0).pt() <= 135.0)){ Manager()->FillHisto("ABCD_C2", Ak15jetConstituents.at(0).size());}

  if( (Ak15jetConstituents.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0) ){ Manager()->FillHisto("ABCD_B1", Ak15jetConstituents.at(0).size());}
  if( (Ak15jetConstituents.at(0).size() >= 14.0) && (Ak15jetConstituents.at(0).size() < 21.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0) ){ Manager()->FillHisto("ABCD_D1", Ak15jetConstituents.at(0).size());}
  if( (Ak15jetConstituents.at(0).size() < 14.0) && (Ak4jets.at(0).pt() <= 220.0) && (Ak4jets.at(0).pt() > 135.0) ){ Manager()->FillHisto("ABCD_D2", Ak15jetConstituents.at(0).size());}

  if( (Ak15jetConstituents.at(0).size() >= 21.0) && (Ak4jets.at(0).pt() > 220.0) ){ Manager()->FillHisto("ABCD_B2", Ak15jetConstituents.at(0).size());}
  if( (Ak15jetConstituents.at(0).size() >= 14.0) && (Ak15jetConstituents.at(0).size() < 21.0) && (Ak4jets.at(0).pt() > 220.0) ){ Manager()->FillHisto("ABCD_E1", Ak15jetConstituents.at(0).size());}
  if( (Ak15jetConstituents.at(0).size() < 14.0)  && (Ak4jets.at(0).pt() > 220.0) ){ Manager()->FillHisto("ABCD_E2", Ak15jetConstituents.at(0).size());}

  Manager()->FillHisto("MET_PT", event.rec()->MET().pt());
  Manager()->FillHisto("MET_ETA", event.rec()->MET().eta());
  Manager()->FillHisto("MET_PHI", event.rec()->MET().phi());
  
  Manager()->FillHisto("Sphericity", sphericity_value);

  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void user::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
}


