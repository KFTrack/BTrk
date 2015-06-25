//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPoca.cc,v 1.46 2006/03/25 15:15:55 brownd Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner, largely taken from Art Snyder
//
//------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/Trajectory.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/BaBar/ErrLog.hh"
using namespace CLHEP;

TrkPoca::TrkPoca(const Trajectory& traj1, double f1, 
		 const Trajectory& traj2, double f2, double prec) 
        : TrkPocaBase(f1,f2,prec), _doca(-9999.) 
{
  TrkPocaTraj ptraj1(traj1,f1,false);
  TrkPocaTraj ptraj2(traj2,f2,false);
  minimize(ptraj1,ptraj2);
  _flt1 = ptraj1._flt;
  _flt2 = ptraj2._flt;
  if (status().failure()) return;
  calcDist(traj1,traj2);
}

TrkPoca::TrkPoca(const Trajectory& traj1, double f1, bool rflt1,
		 const Trajectory& traj2, double f2, bool rflt2,
		 double prec) 
  : TrkPocaBase(f1,f2,prec), _doca(-9999.) 
{ 
  TrkPocaTraj ptraj1(traj1,f1,rflt1);
  TrkPocaTraj ptraj2(traj2,f2,rflt2);
  minimize(ptraj1,ptraj2);
  _flt1 = ptraj1._flt;
  _flt2 = ptraj2._flt;
  if (status().failure()) return;
  calcDist(traj1,traj2);
}



TrkPoca::TrkPoca(const Trajectory& traj, double flt, 
                 const HepPoint& pt, double prec) 
  : TrkPocaBase(flt,prec),_doca(-9999.)
{
  TrkPocaTraj ptraj(traj,flt,false);
  minimize(ptraj,pt);
  _flt1 = ptraj._flt;
  if (status().failure()) return;
  _doca = (traj.position(flt1()) - pt).mag();
}

TrkPoca::TrkPoca(const Trajectory& traj, double flt, bool rflt,
                 const HepPoint& pt, double prec) 
  : TrkPocaBase(flt,prec) ,_doca(-9999.)
{
  TrkPocaTraj ptraj(traj,flt,rflt);
  minimize(ptraj,pt);
  _flt1 = ptraj._flt;
  if (status().failure()) return;
  _doca = (traj.position(flt1()) - pt).mag();
}

TrkPoca::TrkPoca() :
  _doca(0.0)
{}

TrkPoca::TrkPoca(const TrkPoca& other) :
  TrkPocaBase(other),_doca(other._doca)
{}

TrkPoca::~TrkPoca()
{}

TrkPoca*
TrkPoca::clone() const {
  return new TrkPoca(*this);
} 

TrkPoca&
TrkPoca::operator = (const TrkPoca& other) {
  if(this != &other){
    TrkPocaBase::operator =(other);
    _doca = other._doca;
  }
  return *this;
}

void 
TrkPoca::calcDist(const Trajectory& traj1, const Trajectory& traj2) 
{
  // Does a final calculation of the distance -- getting the sign right.
  //  In case of (nearly) parallel, returns error and distance.

  // A bunch of unsightly uninitialized variables:
  static Hep3Vector dir1,dir2;
  static HepPoint pos1,pos2;

  traj1.getInfo(flt1(), pos1, dir1);
  traj2.getInfo(flt2(), pos2, dir2);
  Hep3Vector delta = pos2 - pos1;
  Hep3Vector between = dir1.cross( dir2 );  // cross-product

  if (status().success() != 3) { // Not parallel:
    between = between.unit();
    double size = delta.mag();
    _doca = delta.dot(between) > 0.0 ? size : -size;
  } else {  // Parallel:
    // Arbitrary sign convention
    double sign = (dir1.dot(dir2) > 0.) ? -1. : 1.;
    _doca = (delta - delta.dot(dir1) * dir1).mag();
    _doca *= sign;
  }
}

double
TrkPoca::doca() const{
  return _doca;
}
