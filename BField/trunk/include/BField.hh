//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: BField.hh 497 2010-01-14 09:06:53Z stroili $
//
// Description:
//	The BField class contains the magnetic field map.  This class is now a
//      virtual class
//
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Bob Jacobsen
//      Dave Brown, 12-23-96; Added Kalman filter tracking functions, prepared
//                            for a full field map implementation
//      Fri Sep 25 11:57:57 PDT 1998 - A. Snyder
//      Remove bFieldVectDF and bFieldZDF
//      The are replaced by derivative
//      Dave Brown, 12/20/99  Make this a virtual class, remove the
//      obsolete and dangerous static bFieldNom;
//
// Copyright Information:
//	Copyright (C) 1995, 1996
//
//------------------------------------------------------------------------

#ifndef BFIELD_HH
#define BFIELD_HH

#include "BaBar/Constants.hh"
#include "CLHEP/Geometry/HepPoint.h"
class HepPoint;
#include "CLHEP/Vector/ThreeVector.h"
class DifVector;
class DifPoint;
class DifNumber;
class Code;

class BField {


public:

  BField();
  virtual ~BField() = 0;  // make sure this class is never instantiated

// return the momentum implicit for a given point, direction, and curvature
// Note that all the functions below which return momenta assume the Kalman track
// fit convention that curvature be interpreted in terms of the nominal field.
// If a different physical interpretation of curvature is being used, these
// functions will NOT BE CORRECT.  Conversely, any use of the raw BField
// accesors (bFieldVect in particular) to compute the momentum of a Kalman fit
// given its direction and curvature will NOT BE CORRECT.
  double momentum
  (const HepPoint &position,
   const Hep3Vector &direction,
   double curvature) const;
  
  Hep3Vector momentumVector(const HepPoint&,const Hep3Vector&,
			  double curvature) const;
  // Differential number version
  DifVector momentumDfVector(const HepPoint&,const DifVector&,
			  DifNumber curvature) const;
  //  Invert the above
  double curvature(const HepPoint&,const Hep3Vector& momentum,
			   const double& charge) const;

  // return the magnetic field vector at a point
  virtual Hep3Vector bFieldVect
  (const HepPoint& point = HepPoint(0,0,0)) const = 0;

  //check if a point is ok 
  virtual Code pointOk(const HepPoint& point)const;
 
  //return divergence of field at a point
  double divergence
  (const HepPoint &point = HepPoint(0,0,0),double step=0.01)const;

  //return curl of field at a point
  Hep3Vector curl
  (const HepPoint &point = HepPoint(0,0,0),double step=0.01)const;

  //derivative of field along direction of del
  Hep3Vector derivative
  (const HepPoint& point,
   Hep3Vector& del)const;

  // return the Z component of the magnetic field at a point
  double bFieldZ(const HepPoint& point = HepPoint(0,0,0)) const;

  // return nominal field
  virtual double bFieldNominal()const;

  //print 
  virtual void print(std::ostream& o)const;

  // for converting radius*field to GeV/c of transverse momentum
  static const double cmTeslaToGeVc;

private:
// momentum can be defined without knowledge of the position since the
// Kalman correction implies curvature be interpreted in terms of the
// nominal field only
  double momentum
  (const Hep3Vector& direction,
   double curvature) const;

  bool _cacheHot;
  double _nom;
};

#endif

