//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkDchHotSelector.hh,v 1.1 2001/07/18 00:42:42 brownd Exp $
//
// Description:
//   class TrkDchHotSelector.  Trivial implementation of TrkHotSelector for
//   selecting Dch hots.
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

#ifndef TRKDCHHOTSELECTOR_HH
#define TRKDCHHOTSELECTOR_HH

#include "TrkBase/TrkHotSelector.hh"

class TrkDchHotSelector : public TrkHotSelector{
public:
// only a default constructor
  TrkDchHotSelector(bool activeOnly = true);
  virtual ~TrkDchHotSelector();
// The Function
  virtual bool useHot(const TrkHitOnTrk& hot) const;
private:
  bool _ignoreActive;
// disallow
  TrkDchHotSelector(const TrkDchHotSelector&);
  TrkDchHotSelector& operator = (const TrkDchHotSelector&);
};
#endif
