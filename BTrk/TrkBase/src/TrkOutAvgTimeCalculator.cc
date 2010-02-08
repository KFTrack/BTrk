//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkOutAvgTimeCalculator.cc,v 1.3 2003/01/21 12:55:10 raven Exp $
//
// Description:
//   See header file
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//      Copyright (C) 2001      UC, San Diego
//
// Author List:
//      Gerhard Raven 7/19/01
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkOutAvgTimeCalculator.hh"
#include "TrkBase/TrkHotSelector.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkHitList.hh"
#include <vector>
#include "ErrLogger/ErrLog.hh"


TrkOutAvgTimeCalculator::TrkOutAvgTimeCalculator(const TrkHotSelector& selector, double maxpull) :
  TrkTimeCalculator(selector),_maxpull(maxpull)
{}

TrkOutAvgTimeCalculator::~TrkOutAvgTimeCalculator()
{}

bool
TrkOutAvgTimeCalculator::trackTime(const TrkRecoTrk& trk,
                                   double& time, double& timeerr,
                                   int& nHotsUsed) const
{
  const TrkHitList *hotlist = trk.hits();
  std::vector<ws> l;
  ws sum;
  for(TrkHitList::hot_iterator i(hotlist->begin());i!=hotlist->end();++i){
    double hottime,hottimeerr;
    if(!useHot(*i) || !i->timeResid(hottime,hottimeerr)) continue;
    hottime*=1e9; hottimeerr*=1e9;
    ws x(hottime,double(1)/(hottimeerr*hottimeerr));
    l.push_back(x);
    sum +=x;
  }
  if (_maxpull>0 && sum.n()>2 && sum.isPhysical()) { // do pull based pruning
    std::vector<ws>::iterator worst=l.end();
    do {
        worst=l.end(); double worstPull = _maxpull;
        for(std::vector<ws>::iterator i=l.begin();i!=l.end();++i) {
              ws u(sum); u-=*i;
              if (!u.isPhysical()) {
                      ErrMsg(warning) << "rounding problem?" << endmsg;
                      continue;
              }
              double p = fabs(u.pull( *i ));
              if (p>worstPull) { worstPull = p; worst = i; }
        }
        if (worst!=l.end()) {
                sum-=*worst; l.erase(worst);
        }
    } while (sum.isPhysical() && sum.n()>2 && worst!=l.end());
  }
  if (sum.isPhysical()) {
          time = sum.mean()*1e-9 + trk.trackT0();;
          timeerr = sum.sigma()*1e-9;
          nHotsUsed = sum.n();
          return true;
  }
  return false;
}
