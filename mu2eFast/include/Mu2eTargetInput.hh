// generate mu2e conversions on targets.  Distributions given by
// beam parameters
#ifndef Mu2eTargetInput_HH
#define Mu2eTargetInput_HH
#include "mu2eFast/Mu2eInput.hh"
#include <TRandom3.h>
#include <TLorentzVector.h>

class PdtEntry;
class PacConfig;
class TRandom;
class TGraph;
class PacDetector;

class Mu2eTargetInput : public Mu2eInput {
public:
  // momentum spectrum description
  enum spectrum{flat=0,file};
  // construct from the detector
  Mu2eTargetInput(PacConfig& config, const PacDetector* det);
  ~Mu2eTargetInput();  
  // override virtual interface
  virtual bool nextEvent(Mu2eEvent& event);
  virtual void rewind();
protected:
  TParticle* create(const TLorentzVector& pos, const TLorentzVector& mom) const;
  void createPosition(TLorentzVector& pos) const;
  void createMomentum(TLorentzVector& mom) const;
  
  void prepareBeam(PacConfig& config);
  
private:
  const PdtEntry* _pdt;
  double _beamxsig;
  double _beamtsig;
  double _beamzlambda;
  double _p_min, _p_max;
  double _cost_min, _cost_max;
  double _mass2;
// timing
  double _lifetime;
  double _bunchtime;
  double _lnscale;
  double _lnsigma;
  double _lntheta;
// spectrum type
  spectrum _stype;
  TGraph* _invintspect;
// target description
  std::vector<double> _diskradii;
  std::vector<double> _halfthickness;
// seed
  int _seed;
// event counters
protected:
  mutable TRandom3 _rng;
  std::vector<double> _diskz;
  unsigned _nevents;
  unsigned _ievt;
};
#endif
