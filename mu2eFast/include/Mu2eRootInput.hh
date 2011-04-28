//
//  Input for reading event from ROOT data.  David Brown, LBNL 7 Jan. 2011
//
#ifndef Mu2eRootInput_HH
#define Mu2eRootInput_HH

#include "mu2eFast/Mu2eInput.hh"
#include <Rtypes.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <vector>
#include <string>
class PacConfig;

class Mu2eRootInput : public Mu2eInput {
public:
// create from a single file; the constructor will open the file, and start reading from
// the specified event.  nevents < 0 means read all the events
  Mu2eRootInput(PacConfig& config);
// close the files on destruction
  virtual ~Mu2eRootInput();
// override virtual interface
  virtual bool nextEvent(Mu2eEvent& event);
  virtual void rewind();
// accessors
  const TFile* file() const { return _file; }
  const TTree* tree() const { return _tree; }
// root tree name
  const std::string& treename() { return _treename; }
  void setTreeName(const std::string& tname) { _treename = tname; }
  unsigned nread() const { return _nread; }
protected:
// root structure
  TFile *_file;
  TTree *_tree;
  TBranch *_bparticles, *_bevtnum, *_bevtwt, *_bnevt, *_bnpar;
// tree name; default value can be over-ridden
  std::string _treename;
// branch data; these are controlled by root
  TClonesArray *_particles;
  Int_t _evtnum;
  Float_t _evtwt;
  UInt_t _nevt, _npar;
// allow scaling and shifting the time
  double _tscale;
  double _toffset;
// automatic rewind until specified # of events are read: infinite loop for nevents<0!
  bool _loop;
// initialization information; starting event number
  unsigned _ifirst;
// event counter
  int _nevents;
  unsigned _nread;
  int _ievt;
  int _numevt;
};
#endif
