//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMiniTrk_001.hh,v 1.4 2003/10/02 23:23:02 brownd Exp $
//
// Description:
//      class Kalminitrk.  This adds fully implemented functions to the TrkKalminitrk
//      interface, which allow trivial subsequent implementation of BaBar
//      persistence.  This is still a virtual class.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 11/14/00
//------------------------------------------------------------------------
#ifndef KALMINITRK_001_HH
#define KALMINITRK_001_HH

#include "TrkBase/TrkHelixData_001.hh"
#include "BaBar/BaBarODMGTypes.h"
#include <vector>
class KalMiniRep;
class TrkRecoTrk;
class ComPackFlatFloat;

class KalMiniTrk_001 {
public:
  KalMiniTrk_001(); // needed for Objectivity
  KalMiniTrk_001(const TrkRecoTrk* trk);
// build from raw data, to support persistenc
  KalMiniTrk_001(d_UShort hypomap,d_UShort t0,d_UShort id);
  virtual ~KalMiniTrk_001();
  KalMiniTrk_001& operator = (const KalMiniTrk_001& other);
// TrkKalMiniTrk_001 interface that can be implemented here
  virtual PdtPid::PidType defaultHypo() const;
// the following will return PdtPid::null if the fit failed
  virtual PdtPid::PidType fitHypo(PdtPid::PidType hypo) const;
// count the # of good and failed fits
  unsigned nFit() const;
  unsigned nFailed() const;
  unsigned nMapped() const;
// return a vector of the real fit PdtPids
  void fits(std::vector<PdtPid::PidType>& hypos) const;
// return the index to a given PdtPid, -1 if there's no such fit
  int index(PdtPid::PidType) const;
  virtual unsigned long trackId() const;
  virtual double usedT0() const;
// access to raw data
  d_UShort hypoMap() const { return _hypomap;}
  d_UShort rawT0() const { return _t0;}
  d_UShort rawId() const { return _id; }
  void setData(d_UShort& hypomap,d_UShort& t0,d_UShort& id) const {
    hypomap = _hypomap; t0 = _t0; id = _id; }
private:
// store some non-persistent-specific data members
  d_UShort _hypomap; // map of hypos.  This includes the default hypo
  d_UShort _t0; // t0 (packed)
  d_UShort _id; // id number
// statics
  static const ComPackFlatFloat _packt0; // pack the T0
  static const unsigned _hypomsk; // bit mask for pid hypos
  static const unsigned _hyposft; // bit mask for pid hypos
  static const unsigned _default; // default hypo flag
  static const unsigned _failed; // failed fit flag

  static unsigned hypoShift(unsigned); // utility function
  PdtPid::PidType hypo(unsigned) const;
};

#endif
