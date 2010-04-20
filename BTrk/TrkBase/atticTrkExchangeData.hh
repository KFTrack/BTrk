//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkExchangeData.hh,v 1.1 2003/07/17 20:26:33 brownd Exp $
//
//  Description:
//  Class TrkExchangeData; a compact representation of a helix.  This
//  describes only the helix parameters, not the errors
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Infomation;
//	Copyright (C) 2003	Lawrence Berkeley Laboratory
//
// Author(s): Dave Brown 07/17/03
//
//------------------------------------------------------------------------

#ifndef TRKEXCHANGEDATA_HH
#define TRKEXCHANGEDATA_HH

#include "BaBar/BaBarODMGTypes.h"
#include "TrkBase/TrkExchangePar.hh"

template <class T> class ComPackBase;
class ComPackExpFloat;
class ComPackSignedExpFloat;
class ComPackFlatFloat;

class TrkExchangeData {
public:
// default constructor for Objectivity
  TrkExchangeData();
// construct from a ExchangePar
  TrkExchangeData(const TrkExchangePar*);
// copy and equivalence are OK
  TrkExchangeData(const TrkExchangeData&);
  TrkExchangeData& operator =(const TrkExchangeData&);
// NO VIRTUAL DESTRUCTOR.  Yes, this is on purpose, to cut down on overhead
  ~TrkExchangeData();
// allow separate unpacking of all parameters;
// the helix function _RETURNS OWNERSHIP_
  TrkExchangePar* exchange() const;
// access to raw data
  const d_UShort& parameters(int index ) const {
    return _params[index]; }
private:
  friend class TrkExchangeDataK;
// array of parameter and diagonal covariance terms.
// these are indexed by the same enum defined in TrkExchangePar.
  d_UShort _params[TrkExchangePar::nParam];
  static const ComPackBase<double>& paramPacker(int);
// statics used for packing
  static const ComPackSignedExpFloat _packd0;
  static const ComPackSignedExpFloat _packz0;
  static const ComPackFlatFloat _packphi0;
  static const ComPackSignedExpFloat _packomega;
  static const ComPackFlatFloat _packlambda;
};

#endif
