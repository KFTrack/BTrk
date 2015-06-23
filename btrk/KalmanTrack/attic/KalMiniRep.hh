//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMiniRep.hh,v 1.39 2007/09/04 23:31:11 brownd Exp $
//
//   Description: class KalMiniRep.  Implementation of TrkRep to store a
//   Kalman fit reconstituted from the 'mini'.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 11/6/00
//------------------------------------------------------------------------

#ifndef KALMINIREP_HH
#define KALMINIREP_HH

#include "TrkBase/TrkRep.hh"
#include "ProbTools/ChisqConsistency.hh"
#include "DnaBase/DnaPtr.hh"
#include <vector>

class KalRep;
class TrkSimpTraj;
class KalTrk;
class KalContext;
class TrkDifPieceTraj;
class TrkHotListFull;
class SvtHitOnTrack;

// Class interface //
class KalMiniRep : public TrkRep {
public:
// enum to describe the state of this mini-rep;  In cache mode, all tracking function
// access will fail.  This small-size, small cpu state is created directly from the
// persistent object.  Hot mode operates just like a normal KalRep, while fit mode
// operates by returning stored fit results, properly extended through material and bfield
// effects as necessary.  Note that a mini-rep can be toggled between hot and fit states,
// provided sufficient information exists.  This allows comparing results, or 'upgrading'
// the rep for interior usage (such as kink finding).
  enum miniState { seed=-1,cache=0,extendedcache,hots};
// construct from a list of HOTs, constraints, and/or fit results.
// An explicit seed trajectory is required, though it can be one of the entries on the
// fit or constraint lists.  Note that the KalMiniRep will _TAKE OWNERSHIP_ of the hots,
// constraints, and fit results and seeds on construction.
  KalMiniRep(TrkRecoTrk* trk,PdtPid::PidType hypo,
	     TrkSimpTraj* seed,
	     const KalContext& context,
	     DnaPtr<TrkHotList> fullhotl,
	     TrkHotList* emptyhotl,
	     const std::vector<TrkSimpTraj*>& constraints,
	     const std::vector<TrkSimpTraj*>& extendedcache,
	     double chiprob);
// construct from the default rep.  The allows the KalRep to be cloned.  
// the constraints and fits must still be explicit provided
  KalMiniRep(KalMiniRep* defrep, PdtPid::PidType hypo,
	     const std::vector<TrkSimpTraj*>& constraints,
	     const std::vector<TrkSimpTraj*>& extendedcache,
	     double chiprob);
// copy constructor(s)
  KalMiniRep(const KalMiniRep&,TrkRecoTrk* trk);
  KalMiniRep(const KalMiniRep&,PdtPid::PidType hypo);
  virtual ~KalMiniRep();
  // clone() used to copy tracks; cloneNewHypo() for new hypos within track
  virtual KalMiniRep* clone(TrkRecoTrk* newTrack) const; // covariant return
  virtual TrkRep* cloneNewHypo(PdtPid::PidType hypo);
//******************************************
// TrkRep interface
//******************************************
  virtual HepPoint                  position(double fltL)         const;
  virtual Hep3Vector                direction(double fltL)        const;
  virtual BbrPointErr               positionErr(double fltL)      const;
  virtual BbrVectorErr              directionErr(double fltL)     const;

  virtual TrkDifTraj&       traj();
  virtual const TrkDifTraj& traj()         const;
  virtual int               nDof()         const;
  virtual int               charge()       const;
  virtual Hep3Vector        momentum(double fltL=0.)      const;
  virtual double            pt(double fltL=0.)            const;
  virtual BbrVectorErr      momentumErr(double fltL=0.)      const;
  virtual TrkExchangePar    helix(double fltL=0.) const;
// TrkFit functions
  virtual HepMatrix         posmomCov(double fltL) const;
  virtual void              getAllCovs(double fltL,
				       HepSymMatrix& xxCov,
				       HepSymMatrix& ppCov,
				       HepMatrix& xpCov) const;
  virtual void              getAllWeights(double fltL,
					  HepVector& pos,
					  HepVector& mom,
					  HepSymMatrix& xxWeight,
					  HepSymMatrix& ppWeight,
					  HepMatrix&    xpWeight) const;
// note that a KalInterface can attach to this rep (as well as a KalMiniInterface)
  virtual const IfdKey&     myKey()           const;
// HOT-manipulation functions.  These simply call down to the underlying KalRep (if any).
  virtual void            addHot(TrkHitOnTrk *theHot);
  virtual void            removeHot(TrkHitOnTrk *theHot);
  virtual void            activateHot(TrkHitOnTrk *theHot);
  virtual void            deactivateHot(TrkHitOnTrk *theHot);
  virtual bool            resid(const TrkHitOnTrk *theHot, 
                                double &residual, double &residErr,
                                bool exclude=false) const;
// TrkFitStatus interface.  Override this so both the KalMiniRep and the
// KalRep get the history messages added to the KalMiniRep
  virtual void addHistory(const TrkErrCode& status,const char* modulename);
// extension.  These call down to the KalRep
  TrkErrCode extendThrough(const TrkVolume&,trkDirection trkDir = trkOut ) const;
  TrkErrCode extendThrough(double fltlen) const;
// fitting
  virtual TrkErrCode fit();
// override access to chisq information.  In the mini-rep, consistency is the
// natural variable
  virtual ChisqConsistency    chisqConsistency() const;
  virtual double chisq() const;
//******************************************
// Printing
//******************************************
  virtual void printAll(std::ostream& ostr) const;
  virtual void print(std::ostream& ostr) const;
//******************************************
// mini-rep specific functions
// Printing
//******************************************
// toggle the state.  Note that cache may not be toggled to (that will return
// an error).  The
// return code should _always_ be checked, as the rep may not be able
// to move to the specified mode (in which case the state remains unchanged).
  TrkErrCode changeState(miniState newstate);
// access the current state
  miniState currentState() const { return _state; }
// simple common utility function
  bool active() const { return _state > cache; }
// return the current KalRep
  const KalRep* kalRep() const { return kalRep(_state); }
  KalRep* kalRep() { return kalRep(_state); }
// access the input parameters
  const TrkSimpTraj* seedTrajectory() const { return _seed.get(); }
  const std::vector<TrkSimpTraj*>& fitTrajectories() const { return _fittrajs;}
  const std::vector<TrkSimpTraj*>& constraintTrajectories() const { return _contrajs;}
  void fitTrajectories(std::vector<TrkSimpTraj*>& fits,const char* stream="Default") const;
  void constraintTrajectories(std::vector<TrkSimpTraj*>& fits) const;
// direct access to the underlying KalReps.  This is for expert-only use, as the state
// might not be defined or definable!
  const KalRep* kalRep(miniState state) const;
  KalRep* kalRep(miniState state);
// override the definition of the valid flightlength.  When using cached
// results, this is restricted
  virtual bool validFlightLength(double fltL,double tolerance=0.0)      const;
  virtual void updateHots();
// override hotlist accessor
  virtual TrkHotList*             hotList();
  virtual const TrkHotList*       hotList() const;
// access to the local trajectory of the mini-rep
  const TrkSimpTraj* localTrajectory(double fltlen,double& loclen) const;
  const TrkDifPieceTraj& pieceTraj() const { return *_ptraj; }
// direct access to a cached fit
  const TrkSimpTraj* fitTraj(double fltlen) const;
// give the states names
  static const char* stateName(miniState);
// allow setting the flag whether refit KalRep should override the
// cached fit and constraint values when re-persisting.
  static void setUseKalRep(bool use) { _useKalRep = use; }
  static bool useKalRep() { return _useKalRep; }
// access to the explicit full and empty hot lists (if they exist)
  const TrkHotList* fullHotList() const { return _fullhotlist.rawPtr(); }
  const TrkHotList* emptyHotList() const { return _emptyhotlist.get(); }
// access to context
  const KalContext& kalContext() const { return _kalcon; }
// set static flag for auto-extend
  static void setAutoExtend(bool);
// we can use dEdx to constrain momentum
  bool svtdEdx() const;
private:

  //For borrowing without creating the transient 
  friend class KalMiniTrkK;
  friend class KalMiniRepK_001;
  const DnaPtr<TrkHotList> fullHotListRef() const { return _fullhotlist; }
  void setSvtdEdx(bool usesvtdedx); 

// KalContext
  const KalContext& _kalcon;
// KalReps which actually implement the interace.  The first (may) holds
// the fit result from hots, the 2nd stores extensions of cached trajectories (if any)
  std::auto_ptr<KalRep> _hotrep;
  std::auto_ptr<KalRep> _xrep;
// explicit references to the full and empty hot lists
  DnaPtr<TrkHotList> _fullhotlist;
  std::auto_ptr<TrkHotList> _emptyhotlist;
// owns a copy of the seed trajectory
  std::auto_ptr<TrkSimpTraj> _seed;
  typedef std::vector<TrkSimpTraj*> TrkSimpTrajs;
// array of trajectories to be used in the trajectory-based representation.
  TrkSimpTrajs _fittrajs;
// array of trajectories to be used as constraints in the hot-based representation
  TrkSimpTrajs _contrajs;
// current state
  miniState _state;
// my own trajectory; can respond in cache state (if fit results exist)
  std::auto_ptr<TrkDifPieceTraj> _ptraj;
// my own consistency
  ChisqConsistency _consistency;
  bool _usesvtdedx;
// statics
  static bool _useKalRep; // whether to use KalRep fits or cached fits
  static bool _autoExtend;  // whether to auto-extend fits in cache mode
// create the KalReps
  TrkErrCode createHotKalRep();
  TrkErrCode createXKalRep();
// auto-extension
  void extendIfNeeded(double fltlen) const;
// utility functions
  void setFoundRange();
  void fillHots(std::vector<const SvtHitOnTrack*>&,const TrkHotList* hots) const;
// disallow
  KalMiniRep& operator= (const KalMiniRep&);
  KalMiniRep(const KalMiniRep&);
};
#endif
