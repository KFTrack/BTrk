//
//  Define the input file structure for reading ROOT data particles, and setup
// a structure to read/rewind etc.  David Brown, LBNL 7 Jan. 2011
//
#include "mu2eFast/Mu2eRootInput.hh"
#include "PacEnv/PacConfig.hh"
#include <iostream>
using std::cerr;
using std::endl;

Mu2eRootInput::Mu2eRootInput(PacConfig& config) :
_file(0), _tree(0), _bparticles(0), _bevtnum(0), _bevtwt(0), _bnevt(0), _bnpar(0),
  _treename("Mu2eParticles"),
  _particles(0), _evtnum(0), _evtwt(0), _nevt(0), _npar(0),
  _nread(0),_ievt(0){
// get config information
  const char* filename = config.getcstr("inputfile","none.root");
  _ifirst = config.getint("firstevent",0);
  _nevents = config.getint("nevents",-1);

  _tscale = config.getdouble("timescale",1.0);
  _toffset = config.getdouble("timeoffset",0.0);

  _loop = config.getbool("loop",false);

  _file = TFile::Open(filename,"READ");
  if(_file != 0){
// establish the tree structure
    _tree  = (TTree*)_file->Get(_treename.c_str());
    if(_tree != 0){
// find the branches
      _bparticles = _tree->GetBranch("Particles");
      _bevtnum = _tree->GetBranch("EvtNr");
      _bevtwt = _tree->GetBranch("EvtW");
      _bnevt = _tree->GetBranch("NEvt");
      _bnpar = _tree->GetBranch("NPar");
      _numevt = _tree->GetEntriesFast();
      
      if(_bparticles != 0 && _bevtnum != 0 && _bevtwt != 0 && _bnevt != 0 && _bnpar != 0){
// link branches to local variables
        _bparticles->SetAddress(&_particles);
        _bevtnum->SetAddress(&_evtnum);
        _bevtwt->SetAddress(&_evtwt);
        _bnevt->SetAddress(&_nevt);
        _bnpar->SetAddress(&_npar);
      } else {
        cerr << "Tree branch structure does not correspond to mu2e input in file " << filename << endl;
      }
    } else {
      cerr << "Tree " << _treename << " not found in file " << filename << endl;
      _numevt = 0;
    }
  }
//initialize
  rewind();
}
Mu2eRootInput::~Mu2eRootInput(){
// close the file  
  if(_file != 0){
    _file->Close();
    delete _file;
  }
// also delete local data from branch (??? why doesn't root do this???)
  delete _particles;
}

bool
Mu2eRootInput::nextEvent(Mu2eEvent& event) {
  bool retval(false);
  clear(event,false);
  if(_nread < _nevents){
// move to the next entry in the tree
    int nbytes = _tree->GetEntry(_ievt++);
    if(nbytes >0 && _particles !=0){
      retval = true;
// set particles
      unsigned npar = _particles->GetEntries();
      event._particles.reserve(npar);
      for( int ipar=0; ipar<npar; ipar++ ) {
// must blindly cast:  ugly!
        TParticle* part = (TParticle*)(*_particles)[ipar];
// scale and shift time
        part->SetProductionVertex(part->Vx(),part->Vy(),part->Vz(),_toffset + _tscale*part->T());
        if(part != 0)event._particles.push_back(part);
      }
// update general event information
      event._evtnum = _evtnum;
      event._evtwt = _evtwt;
      event._nevt = _nevt;
      event._npar = _npar;
    }
// reset if we reach the end
    if(_ievt >= _numevt && _loop){
      rewind();
    }
    _nread++;
  }
  return retval;
}

void
Mu2eRootInput::rewind(){
  std::cout << "Rewinding input file " << _file->GetName() << std::endl;
  _ievt = _ifirst;
}
