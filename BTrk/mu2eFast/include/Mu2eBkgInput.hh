//
//  Input for reading event from ROOT data.  David Brown, LBNL 7 Jan. 2011
//
#ifndef Mu2eBkgInput_HH
#define Mu2eBkgInput_HH

#include "mu2eFast/Mu2eRootInput.hh"
#include <Rtypes.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TRandom3.h>
#include <vector>
#include <string>
class PacConfig;

class Mu2eBkgInput : public Mu2eRootInput {
public:
// create from a files, number of events/file to read per 'physics' event
  Mu2eBkgInput(PacConfig& config);
// close the files on destruction
  virtual ~Mu2eBkgInput();
// override virtual interface
  virtual bool nextEvent(Mu2eEvent& event);
private:
  double _bunchtime; // microbunch time spacing
  double _nbkg; // # of background/bunch: includes efficiency
  double _nspread; // fractional spread in stopped muons, interpreted as +-f
  double _norm; // normalization for lifetime sampling
  bool _bkgtime; // use the bkg signal time or not
  double _ymin, _ymax; // limits for sampling bunch fluctuations
  TRandom3 _rng;
};
#endif
