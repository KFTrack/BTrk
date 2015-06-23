//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalPairRep.hh,v 1.20 2003/01/21 13:08:53 raven Exp $
//
// Description:
//      Special subclass of TrkRep used to fit a pair of correlated tracks
//      to a single set of parameters.  This is intended for use in alignment
//      studies, where a single fit to lepton pairs subject to the beam
//      constraint (vertex and boost) is a powerful tool.  KalPairRep inherits
//      most of its functionality (in particular all the fitting
//      mechanism) from KalRep.  Its constructor takes Kalman fits to the 2
//      branches of the pair, plus the constraint which binds them.  Some
//      normal TrkRep functions don't work, for instance:
//      charge().  This is undefined, as the branches are opposte.  KalPairRep
//      returns 0.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 3/21/98
//      Doug Roberts
//------------------------------------------------------------------------
#ifndef KALPAIRREP_HH
#define KALPAIRREP_HH

#include "KalmanTrack/KalRep.hh"
#include "TrkBase/TrkExchangePar.hh"

class BbrDoubleErr;
class BbrPointErr;
class KalPairConduit;

class KalPairRep : public KalRep {

  friend class KalPairConduit;

public:

  KalPairRep(TrkRecoTrk*, const KalRep* seedRep,
	     KalPairConduit* conduit);

  virtual ~KalPairRep();
  
  // Override base class fit function
  TrkErrCode fit();

  // Override nDof functions.  Pair fit has to determine the flightlength
  // at which to perform the transformation, adds an extra parameter
  virtual int nDof() const {return KalRep::nDof(TrkEnums::bothView) - 1;}
  virtual int nDof(TrkEnums::TrkViewInfo view) const
    {return KalRep::nDof(view) - 1;}
  virtual int nDof(double fltlen, trkDirection tdir) const
    {return KalRep::nDof(fltlen, tdir) - 1;}

  bool converged() const;
  //  int iterations() const;
  BbrPointErr prodPoint() const;
  double prodPointChi2() const;
  double pairFltLen();
  KalPairRep* otherPairRep();

  virtual const IfdKey& myKey() const;

private:
  KalPairRep* clone(TrkRecoTrk*) const {return 0;}
  TrkRep* cloneNewHypo(PdtPid::PidType hypo){return 0;}

  KalPairConduit* _conduit;

  // Functions
  void nullConduit() {_conduit = 0;}
};
#endif
