//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkVertex.cc,v 1.6 2005/11/21 22:28:22 ttanabe Exp $
//
// Description:
// Base class for describing TrkRecoTrk-based verteices
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Wouter Hulsbergen
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkVertex.hh"
#include <utility>
#include <algorithm>

TrkVertex::TrkVertex(const PdtEntry* vtype) :
  _status(TrkErrCode::fail),_vtype(vtype)
{}


TrkVertex::~TrkVertex() {}

TrkVertex::TrkVertex(const TrkVertex& other) :
  _usedtracks(other._usedtracks),
  _info(other._info),
  _position(other._position),
  _p4(other._p4),
  _status(other._status),
  _vtype(other._vtype),
  _chisqcon(other._chisqcon)
{}


TrkVertex&
TrkVertex::operator = (const TrkVertex& other) {
  if(this != &other){
    _usedtracks = other._usedtracks;
    _info = other._info;
    _position = other._position;
    _p4 = other._p4;
    _status = other._status;
    _chisqcon = other._chisqcon;
    _vtype = other._vtype;
  }
  return *this;
}

float
TrkVertex::overlap(const TrkVertex& other) const {
// loop over tracks of other vertex
  trkcontainer otrks = other.trackList();
  if(otrks.size() == 0 || _usedtracks.size() == 0) return 0.0;
  unsigned nfound(0);
  for (trkcontainer::const_iterator iother = otrks.begin();
      iother != otrks.end();
      iother++) {
    if(std::find(_usedtracks.begin(),_usedtracks.end(),*iother) != _usedtracks.end())
      nfound++;
  }
  double denom = std::min(otrks.size(),_usedtracks.size());
  return nfound/denom;
}


void
TrkVertex::setConsistency(double chisq, int nDof) {
  _chisqcon = ChisqConsistency(chisq,nDof);
}

bool
TrkVertex::trackInfo(const TrkRecoTrk* trk,TrkVtxInfo& info) const {
  trkinfocontainer::const_iterator ifnd = _info.find(trk);
  bool retval = ifnd != _info.end();
  if(retval)info = ifnd->second;
  return retval;
}

