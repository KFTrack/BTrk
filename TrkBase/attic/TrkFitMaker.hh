//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFitMaker.hh,v 1.28 2004/04/14 04:11:18 brownd Exp $
//
// Description: Abstract base class for FitMaker classes.  FitMakers have 
//   two functions, and only FitMakers are supposed to perform these 
//   functions: 1) creating tracks and 2) changing the TrkRep + Fitter 
//   combination for a track.  
//   Derived classes may have data members -- e.g. cut values 
//   to be set in created fitters.  
//   Derived classes usually implement something like these functions:
//   TrkRecoTrk* makeTrack(const TrkExchangePar& helix, const 
//  			HepAList<TrkHitOnTrk>* hotList, double chi2=-999.);
//   virtual void changeFit(TrkRecoTrk& theTrack) const = 0;
//   
//   But I've found too many exceptions to make these virtual functions.
//   So what this class really does is give friendship access to TrkRecoTrks, 
//   and permits derived classes to muck about to their hearts' content.
//     
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//
//------------------------------------------------------------------------
#ifndef TRKFITMAKER_HH
#define TRKFITMAKER_HH
#include "PDT/PdtPid.hh"
#include <utility>

class TrkHitOnTrk;
class TrkRep;
class TrkRecoTrk;
class TrkExchangePar;
class TrkContext;
class TrkIdManager;
class BField;
class TrkRepIter;
class KalMiniTrkK;

// Class interface //
class TrkFitMaker {

public:
  virtual ~TrkFitMaker();

  void changeDefault(TrkRecoTrk&, PdtPid::PidType) const;

protected:

  // The following functions provide friendship access to TrkRecoTrk 
  //   (including to ctor).
  // gets current TrkReps from track:
  std::pair<TrkRepIter,TrkRepIter> uniqueReps(const TrkRecoTrk& t) const; // was currentReps
  void setRep(TrkRecoTrk&, TrkRep*) const;
  void repointHypo(TrkRecoTrk&, PdtPid::PidType hypo, PdtPid::PidType fit) const;
  void setFitNumber(TrkRecoTrk&, PdtPid::PidType, int) const;
  std::pair<TrkRepIter,TrkRepIter> allReps(const TrkRecoTrk& t) const; // was repPtrs
  TrkRep* getRep(TrkRecoTrk&, PdtPid::PidType) const;
  void addHypoTo(TrkRecoTrk&, TrkRep*, PdtPid::PidType hypo) const;
  TrkRecoTrk* createTrack(PdtPid::PidType, const TrkContext&, double t0) const;
  TrkRecoTrk* createTrack(PdtPid::PidType, long idnum, double t0) const;
// the following functions try to cover the fact that the
// persistence model is completely incompatible with the tracking design.
// They allow direct manipulation of what should be private or constructor-supplied
// arguments.
  void setIdManager(TrkRecoTrk&, TrkIdManager*) const;
  void setBField(TrkRecoTrk&, const BField*) const;
// allow persistence to call these functions
  friend class KalMiniTrkK;
  friend class KalMiniRX;
  friend class KalFit;
};

#endif
