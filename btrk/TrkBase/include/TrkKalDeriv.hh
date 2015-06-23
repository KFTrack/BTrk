//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkKalDeriv.hh,v 1.6 2007/09/24 21:56:27 gapon Exp $
//
// Description:  Pure abstract class adding derrivative functions needed
//  for kalman tracking.  This shouldn't be inherited directly, but is
//  instead brought in either through TrkSimpTraj or TrkDifPieceTraj.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Dave Brown 6/2/97
//
//------------------------------------------------------------------------

#ifndef TRKKALDERIV_HH
#define TRKKALDERIV_HH

#include "CLHEP/Matrix/Matrix.h"

//  Enum defining deflection directions: theta1 is deflection in the
//  rho-z plane, theta2 is deflection in the phi direction.
enum deflectDirection{ theta1=0,theta2};

class TrkKalDeriv {
public:
//  this is a baseclass, therefore it should have a virtual d'tor...
  virtual ~TrkKalDeriv() {} ;
//  Change in param WRT change in direction 
  virtual HepMatrix derivDeflect(double fltlen,deflectDirection idir) const = 0;
// change in param WRT change in transverse position.  This is another
// effect of multiple scattering, only substantial for thick materials
  virtual HepMatrix derivDisplace(double fltlen,deflectDirection idir) const = 0;
//  Similair function for the parameter change WRT change in fractional momentum
  virtual HepMatrix derivPFract(double fltlen) const = 0;

};

#endif
