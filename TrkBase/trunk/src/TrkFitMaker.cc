//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFitMaker.cc,v 1.27 2003/01/21 12:55:09 raven Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkFitMaker.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkRepIter.hh"

TrkFitMaker::~TrkFitMaker()
{
}

std::pair<TrkRepIter,TrkRepIter>
TrkFitMaker::uniqueReps(const TrkRecoTrk& t) const  // was currentReps
{
        return t.uniqueReps();
}

std::pair<TrkRepIter,TrkRepIter>
TrkFitMaker::allReps(const TrkRecoTrk& t) const  // was repPtrs
{
        return t.allReps();
}

void
TrkFitMaker::setRep(TrkRecoTrk& trk, TrkRep* r) const
{
  trk.setRep(r);
}

TrkRep*
TrkFitMaker::getRep(TrkRecoTrk& t, PdtPid::PidType h) const
{
  return t.getRep(h);
}

TrkRecoTrk*
TrkFitMaker::createTrack(PdtPid::PidType hypo, const TrkContext& tc,
                         double t0) const
{
  return new TrkRecoTrk(hypo, tc, t0);
}

TrkRecoTrk*
TrkFitMaker::createTrack(PdtPid::PidType hypo, long idnum,
                         double t0) const
{
  return new TrkRecoTrk(hypo, idnum, t0);
}

void
TrkFitMaker::changeDefault(TrkRecoTrk& t, PdtPid::PidType h) const 
{
  t.changeDefault(h);
}

void
TrkFitMaker::repointHypo(TrkRecoTrk& t, PdtPid::PidType h, PdtPid::PidType f) const
{
  t.repointHypo(h, f);
}

void
TrkFitMaker::setFitNumber(TrkRecoTrk& t, PdtPid::PidType hypo, int newNum) const 
{
  t.setFitNumber(hypo, newNum);
}

void
TrkFitMaker::addHypoTo(TrkRecoTrk& trk, TrkRep* newRep, PdtPid::PidType hypo) const
{
  trk.addHypoTo(newRep, hypo);
}

void
TrkFitMaker::setIdManager(TrkRecoTrk& trk, TrkIdManager* idMan) const
{
  trk.setIdManager(idMan);
}

void
TrkFitMaker::setBField(TrkRecoTrk& trk, const BField* field) const
{
  trk.setBField(field);
}
