//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: GaussConsistency.hh 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Bob Jacobsen, Ed Iskander (LBNL)
//
//      Alexandre Telnov (Princeton), December 2007: 
//         add method "double logErfC(double x)",
//         fix bugs in setting the underflow status; add logLikelihood  
//
// Copyright Information:
//      LBNL, 1996
//      Princeton University, 2007
//
//------------------------------------------------------------------------
#ifndef GAUSSCONSISTENCY_HH
#define GAUSSCONSISTENCY_HH

//-----------------
// BaBar Headers --
//-----------------
#include "BTrk/ProbTools/Consistency.hh"

class GaussConsistency : public Consistency {

  public:

  // default constructor; sets internal state to noMeasure
  GaussConsistency( double delta, double sigma);
  static double logErfC(double x);

  virtual ~GaussConsistency() {}

protected:
  bool calc();      //  consistency of measurement

  double _delta;
  double _sigma;
};

#endif





