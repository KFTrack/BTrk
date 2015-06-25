//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkVertex.hh,v 1.4 2005/11/14 23:48:12 brownd Exp $
//
// Description:
// Base class for describing TrkRep-based verteices
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Wouter Hulsbergen
//
//------------------------------------------------------------------------

#ifndef __TRKVERTEX_HH__
#define __TRKVERTEX_HH__

#include "BTrk/BbrGeom/BbrPointErr.hh"
#include "BTrk/BbrGeom/BbrLorentzVectorErr.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include <vector>
#include <map>

class TrkRep;

class TrkVtxInfo {
public:
  TrkVtxInfo(const double& doca, const double& fltlen) :
    _doca(doca),_fltlen(fltlen){}
  TrkVtxInfo() : _doca(0.0),_fltlen(0.0) {}
  TrkVtxInfo(const TrkVtxInfo& other): _doca(other._doca),_fltlen(other._fltlen){}
  double doca() const { return _doca; }
  double flightLength() const { return _fltlen; }
private:
  double _doca;
  double _fltlen;
};

class TrkVertex {  
public:
  typedef std::vector<const TrkRep*> trkcontainer;
  typedef std::map<const TrkRep*,TrkVtxInfo> trkinfocontainer;

  TrkVertex();
  TrkVertex(const TrkVertex& other);
  TrkVertex& operator = (const TrkVertex& other);
  bool operator == (const TrkVertex& other) { return _usedtracks==other._usedtracks; }
  float overlap(const TrkVertex& other) const;
  virtual TrkVertex* clone() const=0;
  virtual ~TrkVertex();

  int ntracks() const { return _usedtracks.size() ; }
  const trkcontainer& trackList() const { return _usedtracks; }
  double chisq() const { return _chisqcon.chisqValue() ; }
  int nDof() const { return (int)rint(_chisqcon.nDOF()) ; }
  const ChisqConsistency& consistency() const { return _chisqcon; }
  const BbrPointErr& position() const { return _position ; }
  const BbrLorentzVectorErr& p4() const { return _p4 ; }  
  const TrkErrCode& status() const { return _status ; }
// information about particular tracks in the vertex: returns true if track is really in vtx
  bool trackInfo(const TrkRep* trk,TrkVtxInfo& info) const;
protected:
  trkcontainer _usedtracks ;
  trkinfocontainer _info;
  BbrPointErr _position ;
  BbrLorentzVectorErr _p4 ; 
  TrkErrCode _status ;
// to update chisq in subclasses
  void setConsistency(double chisq,int nDof);
private:
  ChisqConsistency _chisqcon;
} ;

#endif
