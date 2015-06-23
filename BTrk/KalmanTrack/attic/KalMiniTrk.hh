//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMiniTrk.hh,v 1.7 2003/08/07 00:19:15 brownd Exp $
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
#ifndef KALMINITRK_HH
#define KALMINITRK_HH

#include "TrkBase/TrkKalTrk.hh"
#include "TrkBase/TrkFitMaker.hh"
#include "TrkBase/TrkHelixData.hh"
#include "BaBar/BaBarODMGTypes.h"
#include <vector>
class KalMiniRep;
class TrkRecoTrk;
class TrkContext;
class KalContext;
class TrkIdManager;
class ComPackFlatFloat;
class TrkKalComposite;
class TrkFitSummary;
class BField;

// This must inherit from TrkFitMaker to access the TrkRecoTrk constructor.
// Too bad, this ties pure transient code into the persistence model!!!??!?!
class KalMiniTrk : public TrkKalTrk, public TrkFitMaker {
public:
  KalMiniTrk(); // needed for Objectivity
  KalMiniTrk(const TrkRecoTrk* trk);
  virtual ~KalMiniTrk();
  KalMiniTrk& operator = (const KalMiniTrk& other);
// TrkKalMiniTrk interface that can be implemented here
  virtual PdtPid::PidType defaultHypo() const;
  virtual PdtPid::PidType fitHypo(PdtPid::PidType hypo) const;
  virtual unsigned long trackId() const;
  virtual double usedT0() const;
  virtual bool isValid(PdtPid::PidType hypo) const;
  virtual bool isCurrent(PdtPid::PidType hypo) const;
  virtual TrkErrCode fitStatus(PdtPid::PidType hypo) const;
  virtual TrkSimpTraj* seedTrajectory() const;
  virtual unsigned nSvt() const;
  virtual unsigned nDch() const;
  virtual unsigned nFit() const;
  int charge(const BField* field) const { return _seed.charge(field); }

// This class can create a TrkRecoTrk from itself, plus some information
// coming from the composite.  This will later be
// called in the 'transient' method of the actual persistent subclass.
// This method _RETURNS OWNERSHIP_.  environment data (BField)
// as well as global event data (TrkIdManager) must be set after construction.
  TrkRecoTrk* createMiniTrack(const KalContext& kalcon,
                              const TrkSimpTraj& seed,
                              TrkHotList* hotlist,
                              std::vector<TrkFitSummary>& fits) const;
// pack the counters.  This function should be private!
  void setCounters(unsigned nsvt,unsigned ndch,unsigned nfit);
// the following is needed due to the composite pattern.  It is not functional
  TrkRecoTrk* transient() const { return 0; }
private:

  friend class KalMiniTrkK;
  friend class TrkKalMiniCompositeK;

// store some non-persistent-specific data members
  d_UShort _hypomap; // map of hypos.  This includes the default hypo
  d_UShort _fitstat; // fit status for each fit
  d_UShort _t0; // t0 (packed)
  d_UShort _id; // id number
  TrkHelixData _seed; // seed
  d_ULong _counters; // counters into svt, dch, fit lists
// utility functions
  void setStatus(KalMiniRep* kmrep) const;
  void fitResult(std::vector<TrkFitSummary>& fits,
                 PdtPid::PidType hypo,TrkEnums::PackFlag flag,
                 std::vector<TrkSimpTraj*>& trajs,double& fitprob) const;
// statics
  static const ComPackFlatFloat _packt0; // pack the T0
  static const unsigned _hypomsk; // bit mask for pid hypos
  static const unsigned _hyposft; // bit mask for pid hypos
  static unsigned hypoShift(unsigned); // utility function
  static const unsigned _validmsk; // bit mask for valid
  static const unsigned _currentmsk; // bit mask for current
  static const unsigned _fitstatmsk; // bit mask for fit status
  static const unsigned _nsvtmsk; // bit masks into counter
  static const unsigned _nsvtsft;
  static const unsigned _ndchmsk;
  static const unsigned _ndchsft;
  static const unsigned _nfitmsk;
  static const unsigned _nfitsft;
};

#endif
