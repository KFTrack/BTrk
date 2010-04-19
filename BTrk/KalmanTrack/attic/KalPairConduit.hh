//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalPairConduit.hh,v 1.5 2001/10/17 14:55:17 steinke Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Doug Roberts
//
//------------------------------------------------------------------------

#ifndef KALPAIRCONDUIT_HH
#define KALPAIRCONDUIT_HH

#include "BbrGeom/BbrDoubleErr.hh"
#include "BbrGeom/BbrVectorErr.hh"
#include "KalmanTrack/KalPairRep.hh"
#include "KalmanTrack/KalPairSite.hh"
#include "KalmanTrack/KalParams.hh"

// Class interface //
class KalPairConduit {

  friend class KalPairRep;
  friend class KalPairSite;

public:
  enum repId {plusRep=0, minusRep};

  KalPairConduit(const TrkDifPieceTraj& plusTraj ,
		 const TrkDifPieceTraj& minusTraj,
		 const BbrVectorErr& beamMom);

  ~KalPairConduit();

  //  bool updateConstraints();

  bool converged() const {return _converged;}
  int iterations() const {return _niter;}

  BbrPointErr prodPoint() const {return _prodPoint;}
  double prodPointChi2() const {return _prodPointChi2;}


protected:	

private:
  // Keep links to KalPairReps and KalPairSites
  KalPairRep* _reps[2];
  KalPairSite* _sites[2];
  int _nreps;
  int _nsites;
  int _niter;

  KalPairRep* otherPairRep(KalPairRep* rep);

  BbrDoubleErr _plusFltLen;
  BbrDoubleErr _minusFltLen;
  BbrPointErr _prodPoint;
  double _prodPointChi2;
  bool _converged;
  bool _needsUpdate;

  // Beam momentum used in the constraint
  BbrVectorErr _beamMom;

  // Actual constraints to be used by KalPairSites
  KalParams* _constraintPar[2];
  KalParams* _inwardPar[2];
  double _fltlens[2];
  bool _siteLAM[2];

  // Functions
  void calcProdPoint();
  bool updateConstraints();

  TrkErrCode coordinateFit();
  void killedBy(KalPairRep* thisRep);

  int thisRepIndex(KalPairRep* thisRep);
  int otherRepIndex(KalPairRep* thisRep);
  int thisSiteIndex(KalPairSite* site);

  void addRep(KalPairRep* newRep);
  void addSite(KalPairSite* newSite, KalPairRep* rep);
  double getFltLen(KalPairRep* thisRep);
  double getFltLen(KalPairSite* thisSite);
  double getFltChi2(KalPairRep* thisRep);
  double getFltChi2(KalPairSite* thisSite);
  KalParams getConstraint(KalPairRep* thisRep);
  KalParams getConstraint(KalPairSite* thisSite);

  bool checkLAM(KalPairSite* site);
  void uploadParams(const KalParams& params, KalPairSite* site);

  // Preempt 
  KalPairConduit&   operator= (const KalPairConduit&);
  KalPairConduit(const KalPairConduit &);
};

#endif







