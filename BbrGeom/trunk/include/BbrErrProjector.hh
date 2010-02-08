//$Id: BbrErrProjector.hh 491 2010-01-13 16:59:16Z stroili $
//Functions to project errors on to planes
//A. Snyder
//SLAC Fri Jul 16 14:51:01 PDT 1999


#ifndef BbrErrProjector_HH
#define BbrErrProjector_HH

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"

class BbrErrProjector {
  
public:

  enum Code {
    ok=0,
    trackParallel,
    uvTooClose,
    matrixWrongSize
  };

  BbrErrProjector() {}

  //project symmetric uncorrelated transverse error on to a plane
  static Code localError
  (double sigma,		// transverse error on track at intersection
   const Hep3Vector& direction,	// track direction
   const Hep3Vector& uHat,	// direction of local coordinate u
   const Hep3Vector& vHat,	// direction of local coordinate v
   HepSymMatrix& error);	// error matrix projected on u-v plane

//private:
  virtual ~BbrErrProjector() {}	// Hide destructor to prevent instantiation
};

#endif
