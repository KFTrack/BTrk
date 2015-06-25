//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHotSelector.hh,v 1.1 2001/07/18 00:42:43 brownd Exp $
//
// Description:
//   class TrkHotSelector.  An abstract base class for deciding whether a
//   particular HOT is useful for a given application.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2001	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 7/17/01
//------------------------------------------------------------------------

#ifndef TRKHOTSELECTOR_HH
#define TRKHOTSELECTOR_HH

class TrkHitOnTrk;

class TrkHotSelector {
public:
// only a default constructor
  TrkHotSelector(){};
  virtual ~TrkHotSelector(){};
// The Function
  virtual bool useHot(const TrkHitOnTrk& hot) const = 0;
private:
// disallow
  TrkHotSelector(const TrkHotSelector&);
  TrkHotSelector& operator = (const TrkHotSelector&);
};
#endif
