// generate muon decay daughters starting with beam muons, causing them to start in the target

#ifndef Mu2eBeamInput_HH
#define Mu2eBeamInput_HH
#include "mu2eFast/Mu2eTargetInput.hh"
#include "mu2eFast/Mu2eRootInput.hh"
#include <TRandom3.h>
#include <TLorentzVector.h>

class PacSimulate;
class PacSimTrack;

class Mu2eBeamInput : public Mu2eTargetInput {
public:
  // construct from configuration and the simulator
  Mu2eBeamInput(PacConfig& config, const PacSimulate* sim);
  ~Mu2eBeamInput();
  // override virtual interface
  virtual bool nextEvent(Mu2eEvent& event);
  virtual void rewind();
protected:
  bool stopsInTarget(const PacSimTrack* strk) const;
private:
  // use the root file input to describe the beam
  Mu2eRootInput* _rinput;
  const PacSimulate* _sim;
  double _lifetime;
  double _bunchtime;
};
#endif
