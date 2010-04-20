//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: PoissonConsistency.hh 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Bob Jacobsen, Ed Iskander
//
// Copyright Information:
//	Copyright (C) 1996
//
//------------------------------------------------------------------------

#ifndef POISSONCONSISTENCY_HH
#define POISSONCONSISTENCY_HH

//-----------------
// BaBar Headers --
//-----------------
#include "ProbTools/Consistency.hh"

class PoissonConsistency : public Consistency {

  public:

  // default constructor; sets internal state to noMeasure
  PoissonConsistency( int nObs, double mu);

  virtual ~PoissonConsistency() {}

private:
  // helper functions

  void calc();      //  cache consistency of measurement

  double getLikelihood( int n, double mu, double logMu ) const; // stash the formula here
  double getLogLikelihood( int n, double mu, double logMu ) const; // stash the formula here
  void   calcGauss() ;
  void   calcDirectSum() ;
  void   calcTails( double logMu ) ;

protected:
  int    _nObs;
  double _mu;
};

#endif





