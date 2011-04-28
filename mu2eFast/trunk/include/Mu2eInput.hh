//
// Base class for generating mu2e FastSim events.  David Brown, LBNL 7 Jan. 2011
//
#ifndef Mu2eInput_HH
#define Mu2eInput_HH

#include <TParticle.h>
#include <vector>
class PacSimTrack;

// simple structure to keep a single event's info coordinated
struct Mu2eEvent {
  std::vector<TParticle*> _particles;
  std::vector<PacSimTrack*> _strks;
  Int_t _evtnum;
  Float_t _evtwt;
  UInt_t _nevt, _npar;
  void append(Mu2eEvent& other){
    _particles.insert(_particles.end(),other._particles.begin(),other._particles.end());
    _strks.insert(_strks.end(),other._strks.begin(),other._strks.end());
    _npar += other._npar;
  }
};

class Mu2eInput {
public:
  Mu2eInput();
  virtual ~Mu2eInput() = 0;
// create/read the next event, and put the relevant particles into a vector.  Note that
// the returned TParticle objects are OWNED BY THIS CLASS and should NOT be deleted or shared
// outisde this event !!!!!
// return value indicates success (end of file, ... as configured)
  virtual bool nextEvent(Mu2eEvent& event) = 0;
// rewind (return to initial conditions)
  virtual void rewind() = 0;
// clear the event, optionally deleting the content
protected:
  void clear(Mu2eEvent& event,bool del=false);
};

#endif
