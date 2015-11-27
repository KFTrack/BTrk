// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalSmear.hh,v 1.5 2007/09/24 18:21:59 gapon Exp $
//
//  Description:
//  Class to describe some loss of information (smearing).
//
// Copyright Information:
//	Copyright (C) 2002	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 2/5/02
//------------------------------------------------------------------------------
#ifndef KALSMEAR_HH
#define KALSMEAR_HH
#include "BTrk/KalmanTrack/KalSite.hh"
#include "CLHEP/Matrix/SymMatrix.h"
//
//  Define the class
//
class KalSmear : public KalSite {
public:
//
//  Constructors; note that the intersection trajectory MUST be a TrkDifPieceTraj
//
  KalSmear(const CLHEP::HepSymMatrix& smear,const TrkDifPieceTraj* ptraj,
	   double glen);
// construct with a smear factor; this will take the current paremeters from the traj 
// and smear them
  KalSmear(const TrkDifPieceTraj* ptraj,
	   double glen,
	   double smearfac);
//  destructor
  virtual ~KalSmear();
//  Fit function
  bool process(const KalSite*,trkDirection);
// update
  virtual bool update(const TrkDifPieceTraj*,double momentum);
//
//  Access
//
  void printAll(std::ostream& os = std::cout) const;
private:
  KalParams _transport; // parameter transport; only the covariance smearing is used
};
#endif
