//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: BField.cc 497 2010-01-14 09:06:53Z stroili $
//
// Description:
//	Class BField; encapsulates the magnetic field.
//
//      See header for more info.  
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Bob Jacobsen		Original Author
//
// Copyright Information:
//	Copyright (C) 1995	Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------

#include <iostream>
#include <fstream>

#include "BField/BField.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Geometry/HepPoint.h"
#include "difAlgebra/DifNumber.hh"
#include "difAlgebra/DifVector.hh"
#include "difAlgebra/DifPoint.hh"
#include "difAlgebra/Code.hh"
using std::cout;
using std::endl;
using std::fstream;
using std::ostream;

// RKK: this is now badly named.
// Units are now:  mm, ns, MeV, T
const double BField::mmTeslaToMeVc = Constants::c/1.0E3;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//-----------------------------------------------------------------------
// dummy constructor
BField::BField()
        :_cacheHot(false),_nom(0)
{}

// destroy
BField::~BField() {}	    



double
BField::bFieldZ(const HepPoint& p) const
{
  return bFieldVect(p).z();
}

//-----------------------------------------------------------------------
//  Compute the momentum implicit for a trajectory at a specified flight length
//  This is a trivial implementation for a constant field, a full field map
//  would need to inquire of the track where it is in the volume to get the
//  right field magnitude and direction; However, in the case of the Kalman filter
//  even in the presence of a non-trivial field this trivial implementation is 
//  still correct: all of the non-trivial work is done in the Kalman filter...
double
BField::momentum(const HepPoint& where, const Hep3Vector& direction, 
                 double curvature) const 
{
  return momentum(direction,curvature);
}

double
BField::momentum(const Hep3Vector& direction, double curvature) const 
{
  double cosdip=direction.perp();
  return  fabs(mmTeslaToMeVc*bFieldNominal()*cosdip/curvature);
}

Hep3Vector
BField::momentumVector(const HepPoint& point,const Hep3Vector& direction,
                       double curvature) const 
{
  Hep3Vector mom = direction;
  mom.setMag(momentum(direction,curvature));
  return mom;
}

DifVector 
BField::momentumDfVector(const HepPoint& where, const DifVector& direction,
			 DifNumber curvature) const {
//-----------------------------------------------------------------------
//  A real field map would look up the magnitude and direction here

//  This assumes the curvature is defined as the magnitude of the 2nd
//  derrivative of the position with respect to 3-d flight length

  DifVector mom = direction;

  //  DifNumber cosdip = direction.perp();  // assumes direction = unit vector

  DifVector field(bFieldVect(where));
  DifNumber sindip=field.unit()*direction;
  DifNumber cosdip=sqrt(1.0-sindip*sindip);

  DifNumber mag = mmTeslaToMeVc*bFieldNominal()*cosdip/curvature;
  mag.absolute();
  mom *= mag;
  return mom;
}

//  Inverse function to the above, with similair caveats on the implementation
//-----------------------------------------------------------------------
double
BField::curvature(const HepPoint& where,const Hep3Vector& momentum,
                  const double& charge) const 
{
  Hep3Vector field=bFieldVect(where);
  double sindip=field.unit()*momentum.unit();
  double cosdip=sqrt(1.0-sindip*sindip);

  return -charge*mmTeslaToMeVc*bFieldNominal()*cosdip/momentum.mag();
}

double 
BField::divergence (const HepPoint &point,double step) const 
{
  Hep3Vector dx(step,0,0);
  Hep3Vector dy(0,step,0);
  Hep3Vector dz(0,0,step);
  double dBxdX=derivative(point,dx).x();
  double dBydY=derivative(point,dy).y();
  double dBzdZ=derivative(point,dz).z();
  return dBxdX+dBydY+dBzdZ;
}

Hep3Vector
BField::curl(const HepPoint &point,double step) const 
{
  Hep3Vector dx(step,0,0);
  Hep3Vector dy(0,step,0);
  Hep3Vector dz(0,0,step);
  double f=1.0/(2.0*step);
  Hep3Vector dBdX=derivative(point,dx);
  Hep3Vector dBdY=derivative(point,dy);
  Hep3Vector dBdZ=derivative(point,dz);
  return 
    Hep3Vector(dBdY.z()-dBdZ.y(),dBdZ.x()-dBdX.z(),dBdX.y()-dBdY.x());
}

Hep3Vector
BField::derivative (const HepPoint& point, Hep3Vector& del)const
{
  long loop(1);
  while(loop) {
    if(pointOk(point+del).fail()) {
      del=del*0.5;
      continue;
    }
    else if(pointOk(point-del).fail()) {
      del=del*0.5;
      continue;
    }
    loop++;
    if(loop<7) continue;	// too close to edge?
    del=Hep3Vector(0,0,0);
    return Hep3Vector(0,0,0);
  }
  double f=2.0*del.mag();
  return (bFieldVect(point+del)-bFieldVect(point-del))*f;
}

void BField::print(ostream &o)const
{
  cout << "Default constant field" << endl;
  cout << "bFieldVect(origin):" << bFieldVect(HepPoint(0,0,0)) << endl; 
  cout << "nominal Field: " << bFieldNominal() << endl;
}

double
BField::bFieldNominal()const
{
  if (!_cacheHot) {
          BField *self=const_cast<BField *>(this);
          self->_nom=bFieldZ();
          self->_cacheHot=true;
  }
  return _nom;
}

Code
BField::pointOk(const HepPoint& aPoint)const
{
  return Code(1,0);
}

Hep3Vector
BField::bFieldVect(const HepPoint& point) const 
{
  int YOU_SHOULD_NEVER_CALL_THIS_FUNCTION_DIRECTLY(-1);
  assert (YOU_SHOULD_NEVER_CALL_THIS_FUNCTION_DIRECTLY == 0);
  return Hep3Vector(0.0,0.0,0.0);
}
