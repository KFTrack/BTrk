// generate mu2e conversions on targets.  Distributions given by
// beam parameters
#ifndef Mu2eTargetInput_HH
#define Mu2eTargetInput_HH
#include "mu2eFast/Mu2eInput.hh"
#include <TRandom3.h>

class PdtEntry;
class PacConfig;
class TRandom;
class GTrack;
class GVertex;

class Mu2eTargetInput : public Mu2eInput {
public:
  // construct from configuration object
  Mu2eTargetInput(PacConfig& config);
  // override virtual interface
  virtual bool nextEvent(Mu2eEvent& event);
  virtual void rewind();
protected:
  TParticle* create();
private:
  const PdtEntry* _pdt;
  double _beamxsig;
  double _beamtsig;
  double _beamzlambda;
  double _pconv;
  double _cost_min, _cost_max;
  unsigned _nevents;
  TRandom3 _rng;
// target description
  std::vector<double> _diskradii;
  std::vector<double> _diskz;
// seed
  int _seed;
// event counters
  unsigned _ievt;
};
#endif
