// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkStoreHypo.hh,v 1.2 2002/10/12 00:07:11 brownd Exp $
//
//  Description: TrkStoreHypo.  Trivial class to record which track
//  fits sampled where should be stored when a TrkRecoTrk is persisted
//
// Copyright Information:
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 11/20/00
//------------------------------------------------------------------------------

#ifndef TRKSTOREHYPO_HH
#define TRKSTOREHYPO_HH

#include "BaBar/PdtPid.hh"

class TrkStoreHypo {
public:
// constructors and destructor
  TrkStoreHypo() : _hypo(PdtPid::null), _fltlen(-9999.0)
  {}
  TrkStoreHypo(PdtPid::PidType hypo,double fltlen) : _hypo(hypo),_fltlen(fltlen)
  {}
  TrkStoreHypo(const TrkStoreHypo& other) : _hypo(other._hypo),_fltlen(other._fltlen)
  {}
  ~TrkStoreHypo(){}
// accessors
  const PdtPid::PidType& hypo() const { return _hypo; }
  double flightLength() const { return _fltlen; }
// rw functions
  TrkStoreHypo& operator = (const TrkStoreHypo& other) {
    if(this != &other){
      _hypo = other._hypo;
      _fltlen = other._fltlen;
    }
    return *this;
  }
  bool operator == (const TrkStoreHypo& other) const;
  bool operator < (const TrkStoreHypo& other) const;
private:
  PdtPid::PidType _hypo;
  double _fltlen;
  static double _flttol;
};

#endif
