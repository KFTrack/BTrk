// generate mu2e conversions on targets.  Distributions given by
// beam parameters
#ifndef Mu2eTargetInput_HH
#define Mu2eTargetInput_HH
#include "mu2eFast/Mu2eInput.hh"
#include <TRandom3.h>

class PdtEntry;
class PacConfig;
class TRandom;
class TGraph;
class GTrack;
class GVertex;
class PacDetector;

class Mu2eTargetInput : public Mu2eInput {
public:
  // momentum spectrum description
  enum spectrum{flat=0,file};
  // construct from configuration object
  Mu2eTargetInput(PacConfig& config);
  // construct from the detector
  Mu2eTargetInput(PacConfig& config, const PacDetector* det);
  // override virtual interface
  virtual bool nextEvent(Mu2eEvent& event);
  virtual void rewind();
protected:
  TParticle* create();
  void prepareBeam(PacConfig& config);
  
private:
  const PdtEntry* _pdt;
  double _beamxsig;
  double _beamtsig;
  double _beamzlambda;
  double _p_min, _p_max;
  double _cost_min, _cost_max;
  unsigned _nevents;
  TRandom3 _rng;
// spectrum type
  spectrum _stype;
  TGraph* _invintspect;
// target description
  std::vector<double> _diskradii;
  std::vector<double> _diskz;
  std::vector<double> _halfthickness;
// seed
  int _seed;
// event counters
  unsigned _ievt;
};
#endif
