#ifndef analysis_cms_exo_23_003_CRDY_h
#define analysis_cms_exo_23_003_CRDY_h
#include <algorithm> // for std::sort
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"
#include "SampleAnalyzer/Interfaces/root/RootMainHeaders.h"

namespace MA5
{
class cms_exo_23_003_CRDY : public AnalyzerBase
{
  INIT_ANALYSIS(cms_exo_23_003_CRDY, "cms_exo_23_003_CRDY")

  public : 
      MAbool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
      void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
      MAbool Execute(SampleFormat& sample, const EventFormat& event);

  private : 
};
}

#endif