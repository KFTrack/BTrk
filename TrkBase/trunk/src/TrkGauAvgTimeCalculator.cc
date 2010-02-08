//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkGauAvgTimeCalculator.cc,v 1.6 2004/09/10 18:00:18 bartoldu Exp $
//
// Description:
//   class TrkGauAvgTimeCalculator.  A simple implementation of TrkTimeCalculator
//   that takes the weighted average of all selected hot's times assuming a
//   Gaussian error.
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

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkGauAvgTimeCalculator.hh"
#include "TrkBase/TrkHotSelector.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkHitList.hh"
#include <math.h>
#include <iterator>

TrkGauAvgTimeCalculator::TrkGauAvgTimeCalculator(const TrkHotSelector& selector) :
  TrkTimeCalculator(selector)
{}

TrkGauAvgTimeCalculator::~TrkGauAvgTimeCalculator(){};


bool
TrkGauAvgTimeCalculator::trackTime(const TrkRecoTrk& trk,
                                   double& time,
                                   double& timeerr,
                                   int& nHotsUsed) const 
{
  double hottime,hottimeerr;
  double wtsum(0.0);
  double wtavg(0.0);
  int n(0);
// get the hot list
  const TrkHitList *hotlist = trk.hits();
  for(TrkHitList::hot_iterator i(hotlist->begin());i!=hotlist->end();++i){
    if(useHot(*i) && i->timeAbsolute(hottime,hottimeerr)){
      double wt = 1.0/(hottimeerr*hottimeerr);
      wtsum += wt;
      wtavg += wt*hottime;
      ++n;
    }
  }
  if(wtsum > 0.0){
    time = wtavg/wtsum;
    timeerr = sqrt(1.0/wtsum);
    nHotsUsed=n;
    return true;
  } else
    return false;
}
