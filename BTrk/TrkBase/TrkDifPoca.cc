//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkDifPoca.cc,v 1.30 2008/03/18 21:43:20 brownd Exp $
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
#include "BTrk/TrkBase/TrkDifPoca.hh"
#include "BTrk/TrkBase/TrkDifTraj.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/difAlgebra/DifVector.hh"
#include "BTrk/difAlgebra/DifNumber.hh"
#include "BTrk/difAlgebra/DifPoint.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/BaBar/ErrLog.hh"
using namespace CLHEP;

TrkDifPoca::TrkDifPoca(const TrkDifTraj& traj1, double f1, 
                       const Trajectory& traj2, double f2, double prec) 
        : TrkPocaBase(f1,f2,prec), _doca(-9999.,0) 
{
  TrkPocaTraj ptraj1(traj1,f1,false);
  TrkPocaTraj ptraj2(traj2,f2,false);
  minimize(ptraj1,ptraj2);
  _flt1 = ptraj1._flt;
  _flt2 = ptraj2._flt;
  if (status().failure()) return;
  calcDist(traj1,traj2);
}


TrkDifPoca::TrkDifPoca(const TrkDifTraj& traj, double flt, 
                       const HepPoint& pt, double prec) 
        : TrkPocaBase(flt,prec), _doca(-9999.,0) 
{ 
  TrkPocaTraj ptraj(traj,flt,false);
  minimize(ptraj,pt);
  _flt1 = ptraj._flt;
  if (status().failure()) return;
  calcDist(traj,pt);
}

TrkDifPoca::TrkDifPoca(const TrkDifPoca& other) :
  TrkPocaBase(other), _doca(other._doca)
{}

TrkDifPoca::~TrkDifPoca()
{}

TrkDifPoca*
TrkDifPoca::clone() const {
  return new TrkDifPoca(*this);
}


void 
TrkDifPoca::calcDist(const TrkDifTraj& traj1, const Trajectory& traj2) 
{
  // Does a final calculation of the distance -- better behaved derivs than 
  //  stepTowardPoca for zero dist.  In case of (nearly) parallel, returns 
  //  distance calculated at whatever point we happen to be at.  
  //  Derivatives are taken (of dist) w/r/t traj2 
  //  parameters, using DifNumbers for the calculation.

  // A bunch of unsightly uninitialized variables:
  static DifVector dir1;
  static DifPoint pos1;
  static Hep3Vector dir2;
  static HepPoint pos2;

  // Request DifNumber info from traj1, ordinary info from traj2
  traj2.getInfo(flt2(), pos2, dir2);
  traj1.getDFInfo2(flt1(), pos1, dir1);

  DifVector delta =   DifPoint(pos2) - pos1;
  Hep3Vector del = pos2-pos1.hepPoint();

  if (status().success() != 3) { // Not parallel:
// here we must go slow but sure
    DifVector between = cross( dir1, dir2 );  // cross-product
    between.normalize();
    _doca = delta * between;
  } else {  // Parallel: Arbitrary sign convention
    _doca = (delta - delta.dot(dir1) * dir1).length();
    if  (dir1.dot(dir2) > 0.) _doca.flipsign();
    if (fabs(_doca.number()) < 0.0001 * precision()) {
      // Parallel and on top of each other (derivatives singular)
      _doca = 0;
      _status.setFailure(3);
    }
  }
}

void 
TrkDifPoca::calcDist(const TrkDifTraj& traj, const HepPoint& pt) 
{
  // Does a final calculation of the distance -- and handles singularity 
  // in derivs when d = 0.

  // A bunch of unsightly uninitialized variables:
  static DifVector dir;
  static DifPoint posTraj;

  DifPoint posPoint(pt);
  traj.getDFInfo2(flt1(), posTraj, dir);  // call faster one, if exists
  DifVector delta = posTraj - posPoint;
  // delta -= dir*delta;  // note: dir*delta is zero, but the derivatives are NOT

  DifNumber dist = delta.length();
  if (dist.number() > 0.01 * precision()) { // d != 0
    _doca = dist;
  } else {
    // d = 0.  Fudge like mad: pick a direction (not parallel to traj) and 
    //   move the point slightly.  This should not happen very often.
    Hep3Vector fudgeVec(0., 0.1 * precision(), 0.0);
    if (dir.x.number() == 0.0 && dir.z.number() == 0.0) {
      fudgeVec = Hep3Vector(0.1 * precision(), 0., 0.);
    }
    posPoint += fudgeVec;
    delta = posTraj - posPoint;
    _doca = dist;
    _status.setSuccess(20, "TrkDifPoca: distance zero, calculation fudged.");
  }
}

double
TrkDifPoca::doca() const {
  return _doca.number();
}
