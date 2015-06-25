//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: AsymGaussConsistency.hh 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//      Consistency subclass for an asymmetric gaussian distribution
//      (different sigmas on two sides of the peak)
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Yury Kolomensky
//      Alexandre Telnov (Princeton), December 2007: 
//         fix bugs in setting the underflow status; add logLikelihood  
//
// Copyright Information:
//	Copyright (C) 1998       Caltech
//	Copyright (C) 2007       Princeton
//
//------------------------------------------------------------------------
#ifndef ASYMGAUSSCONSISTENCY_HH
#define ASYMGAUSSCONSISTENCY_HH
#include "BTrk/BaBar/BaBar.hh"

//-----------------
// BaBar Headers --
//-----------------
#include "BTrk/ProbTools/Consistency.hh"

class AsymGaussConsistency : public Consistency {
public:
  AsymGaussConsistency( double delta, double sigmaMinus, double sigmaPlus);

  virtual ~AsymGaussConsistency() {}

protected:
  bool calc();      //  consistency of measurement

  double _delta;
  double _sigmaMinus;
  double _sigmaPlus;
};

#endif





