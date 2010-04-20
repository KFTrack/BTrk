//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkSvtHotSelector.hh,v 1.1 2001/07/18 00:42:43 brownd Exp $
//
// Description:
//   class TrkSvtHotSelector.  Trivial implementation of TrkHotSelector for
//   selecting Svt hots.
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

#ifndef TRKSVTHOTSELECTOR_HH
#define TRKSVTHOTSELECTOR_HH

#include "TrkBase/TrkHotSelector.hh"

class TrkSvtHotSelector : public TrkHotSelector{
public:
// only a default constructor
  TrkSvtHotSelector(bool activeOnly = true);
  virtual ~TrkSvtHotSelector();
// The Function
  virtual bool useHot(const TrkHitOnTrk& hot) const;
private:
  bool _ignoreActive;
// disallow
  TrkSvtHotSelector(const TrkSvtHotSelector&);
  TrkSvtHotSelector& operator = (const TrkSvtHotSelector&);
};
#endif
