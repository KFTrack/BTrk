// simple single-particle generator for testing
#ifndef Mu2eSimpleInput_HH
#define Mu2eSimpleInput_HH
#include "mu2eFast/Mu2eInput.hh"
#include <TRandom3.h>

class PdtEntry;
class PacConfig;
class TRandom;
class GTrack;
class GVertex;

class Mu2eSimpleInput : public Mu2eInput {
public:
  // construct from configuration object
  Mu2eSimpleInput(PacConfig& config);
  // override virtual interface
  virtual bool nextEvent(Mu2eEvent& event);
  virtual void rewind();
protected:
  TParticle* create();
private:
  const PdtEntry* _pdt;
  double _p_min;
  double _p_max;
  double _cost_min;
  double _cost_max;
  double _r0_mean;
  double _r0_sigma;
  double _z0_mean;
  double _z0_sigma;
  unsigned _nevents;
  TRandom3 _rng;
// seed
  int _seed;
// event counters
  unsigned _ievt;
};
#endif
