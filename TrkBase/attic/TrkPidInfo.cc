//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPidInfo.cc 64 2009-03-27 09:00:51Z stroili $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkPidInfo.hh"
#include "PDT/PdtEntry.hh"
#include "AbsPid/PidMassHypo.hh"
#include "ProbTools/Consistency.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include <iostream>
using std::endl;
using std::ostream;

//------------------------------------------------------------------------
TrkPidInfo::TrkPidInfo(const TrkRecoTrk* trk) : _trk(trk) {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
TrkPidInfo::TrkPidInfo(const TrkPidInfo &rhs) : AbsPidInfo(rhs) , _trk(rhs._trk) {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
TrkPidInfo::~TrkPidInfo() {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
TrkPidInfo&   
TrkPidInfo::operator= (const TrkPidInfo& rhs) {
//------------------------------------------------------------------------
  AbsPidInfo::operator=(rhs);
  _trk = rhs._trk;
  return *this;
}

//------------------------------------------------------------------------
TrkPidInfo* 
TrkPidInfo::clone() const {
//------------------------------------------------------------------------
  return new TrkPidInfo(*this);
}

//------------------------------------------------------------------------
const TrkRecoTrk*
TrkPidInfo::recoTrack() const {
//------------------------------------------------------------------------
  return _trk;
}

//------------------------------------------------------------------------
void
TrkPidInfo::print(ostream& ostr) const {
//------------------------------------------------------------------------
  ostr << "----------------------------------------------------------" << endl;
  ostr << "TrkPid info for trk " << recoTrack()->id() << endl;
  
  int index = 1;
  while ( massHypo(index) != 0) {
    const PidMassHypo* mhyp = massHypo(index);
    ostr << "    Pdt entry: " << mhyp->pdtEntry()->name() << endl;
    ostr << "    ";
    mhyp->consistency().print(ostr);

    index++;
  }
  ostr << "----------------------------------------------------------" << endl;
}

