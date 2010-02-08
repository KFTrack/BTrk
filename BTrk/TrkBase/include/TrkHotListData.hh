//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHotListData.hh,v 1.6 2004/03/24 01:54:43 brownd Exp $
//
//  Description:
//  Class TrkHotListData; a compact representation of a Hot List (without the hots)
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Infomation;
//	Copyright (C) 2002	Lawrence Berkeley Laboratory
//
// Author(s): Dave Brown 10/202/02
//
//------------------------------------------------------------------------

#ifndef TRKHOTLISTDATA_HH
#define TRKHOTLISTDATA_HH

#include "BaBar/BaBarODMGTypes.h"

template <class T> class ComPackBase;
class ComPackExpFloat;
class ComPackSignedExpFloat;
class ComPackFlatFloat;
class ComPackInt;
class TrkHotList;
class TrkHotListEmpty;
class TrkView;

class TrkHotListData {
public:
// default constructor for Objectivity
  TrkHotListData();
// construct from a TrkHotList
  TrkHotListData(const TrkHotList&);
// construct from raw data; this is used for presistence reconstitution
  TrkHotListData(d_UShort svtdata,d_UShort dchdata,d_ULong rangedata);
// copy and equivalence are OK
  TrkHotListData(const TrkHotListData&);
  TrkHotListData& operator =(const TrkHotListData&);
// NO VIRTUAL DESTRUCTOR.  Yes, this is on purpose, to cut down on overhead
  ~TrkHotListData();
// allow separate unpacking of all parameters;
// the hotlist function _RETURNS OWNERSHIP_
  TrkHotListEmpty* hotList() const;
// access the raw data
  d_UShort svtData() const { return _svtdata; }
  d_UShort dchData() const { return _dchdata; }
  d_ULong rangeData() const { return _fndrng; }
  void setData(d_UShort& svtdata,d_UShort& dchdata,d_ULong& rangedata) const {
    svtdata = _svtdata; dchdata = _dchdata; rangedata = _fndrng; }
// processed accessors

  unsigned nPhi() const;
  unsigned nZ() const;
  unsigned nAxial() const;
  unsigned nStereo() const;
  unsigned firstDch() const;
  unsigned lastDch() const;
  TrkView svtView(int layer) const;
  double foundStart() const;
  double foundRange() const;

private:
// Packed data
  d_UShort _svtdata; //# and pattern of svt hots
  d_UShort _dchdata; //# of dch hots
  d_ULong _fndrng; // flight range + dch layer range
// statics used for packing
  static const ComPackSignedExpFloat _packfndstart;
  static const ComPackFlatFloat _packfndrange;
// svt counting + pattern shifts+masks
  static const unsigned _vsft;
  static const unsigned _vmsk;
  static const unsigned _nphisft;
  static const unsigned _nphimsk;
  static const unsigned _nzsft;
  static const unsigned _nzmsk;
// dch counting shifts+masks
  static const unsigned _naxsft;
  static const unsigned _naxmsk;
  static const unsigned _nstsft;
  static const unsigned _nstmsk;
// flight range shifts+masks
  static const unsigned _startsft;
  static const unsigned _startmsk;
  static const unsigned _rangesft;
  static const unsigned _rangemsk;
  static const unsigned _flaysft;
  static const unsigned _flaymsk;
  static const unsigned _llaysft;
  static const unsigned _llaymsk;
// simple counting
  static const unsigned _nsvt;
  static const unsigned _ndch;
};

#endif
