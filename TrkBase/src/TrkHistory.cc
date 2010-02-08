//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHistory.cc,v 1.2 2004/08/06 06:31:42 bartoldu Exp $
//
// Description:  TrkHistory; simple class to record the history
// of a track fit
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Dave Brown
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkHistory.hh"
using std::ostream;

TrkHistory::TrkHistory() : _mod("Unknown")
{}

TrkHistory::TrkHistory(const TrkErrCode& err,const char* modname):
  _err(err),_mod(modname)
{}

TrkHistory::TrkHistory(const TrkHistory& other) : _err(other._err),
						  _mod(other._mod)
{}

TrkHistory&
TrkHistory::operator =(const TrkHistory& other){
  if(&other != this){
    _err = other._err;
    _mod = other._mod;
  }
  return *this;
}

TrkHistory::~TrkHistory()
{}

bool
TrkHistory::operator == (const TrkHistory& other) const {
  return false;
}

void
TrkHistory::print(ostream& os) const {
  os << module() << ": " << status();
}

ostream& operator<<(ostream& os, const TrkHistory& history) {
  history.print(os);
  return os;
}

