//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchHotSelector.hh 89 2010-01-14 12:34:14Z stroili $
//
// Description:
//   class DchHotSelector.  implementation of TrkHotSelector for
//   selecting Dch hots.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2001	UC, San Diego
//
// Author List:
//      Gerhard Raven 7/19/01
//------------------------------------------------------------------------

#ifndef DCHHOTSELECTOR_HH
#define DCHHOTSELECTOR_HH

#include "TrkBase/TrkHotSelector.hh"

class DchHotSelector : public TrkHotSelector{
public:
// only a default constructor
  DchHotSelector(bool activeOnly = true,
                 double minDoca =0, double maxDoca =0);
  virtual ~DchHotSelector();
// The Function
  virtual bool useHot(const TrkHitOnTrk& hot) const;
private:
  bool _ignoreActive;
  double _minDoca ,_maxDoca ;
// disallow
  DchHotSelector(const DchHotSelector&);
  DchHotSelector& operator = (const DchHotSelector&);
};
#endif
