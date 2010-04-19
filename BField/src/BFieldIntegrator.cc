//--------------------------------------------------------------------------
//$Id: BFieldIntegrator.cc 497 2010-01-14 09:06:53Z stroili $
//Description:Class Implementation for |BFieldIntegrator|
//Author List:A. Snyder
//Copyright (C) 1998	SLAC
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "BField/BField.hh"
#include "BbrGeom/Trajectory.hh"
#include "BField/BFieldIntegrator.hh"
#include <vector>

BFieldIntegrator::BFieldIntegrator
(const BField &f):
  _field(f),
  _bNominal(0.0,0.0,f.bFieldNominal()),
  _tolerance(0.1),
  _stepFrac(0.1),_stepCeiling(5.0),_pathMin(0.5),
  _divTolerance(0.1),
  _divStepFrac(0.1),_divStepCeiling(5.0),_divPathMin(0.5)
{
  //  cerr << "BFieldIntegrator: _bNominal = " << _bNominal << endl;
  //  f.print(cerr);
}

BFieldIntegrator::~BFieldIntegrator(){}

Hep3Vector
BFieldIntegrator::deltaMomentum
(const Trajectory *traj, double range[2])const{
  double slo=range[0];
  double sup=range[1];
  return deltaMomentum(traj,slo,sup);
}

Hep3Vector
BFieldIntegrator::deltaMomentum
 (const Trajectory *traj, double slo, double sup)const{

  HepPoint beg=traj->position(slo);
  HepPoint end=traj->position(sup);
  double deltas=sup-slo;
  double smid=(slo+sup)/2;
  double test=deltas*(_field.bFieldVect(end)-_field.bFieldVect(beg)).mag();
  if(fabs(test)<=_tolerance
     ||
     deltas<pathMin() )
    {
      HepPoint mid=traj->position(smid);
      Hep3Vector bMid=field().bFieldVect(mid);
      bMid-=_bNominal;
      Hep3Vector dir=traj->direction(smid);
      Hep3Vector temp=dir*deltas;
      temp*=BField::cmTeslaToGeVc;
      return temp.cross(bMid);
    }else {
      return 
 	deltaMomentum(traj,slo,smid)
 	+deltaMomentum(traj,smid,sup);
    }
}

void 
BFieldIntegrator::divideRange
(const Trajectory *traj,
 double range[2],
 std::vector<double>& posList) const
{
  double slo=range[0];
  double sup=range[1];
  divideRange(traj,slo,sup,posList);
}

void 
BFieldIntegrator::divideRange
(const Trajectory *traj,
 double slo,
 double sup,
 std::vector<double>& posList)const
{
  // trivial test first
  double deltas=sup-slo;
  if(deltas > divPathMin()) {
    double smid=(slo+sup)/2.;
    HepPoint mid=traj->position(smid);
    Hep3Vector bmid(_field.bFieldVect(mid));
    //    HepPoint beg=traj->position(slo);
    //    HepPoint end=traj->position(sup);
    double test = deltas*( bmid - _bNominal).mag();
    //    double change = (_field.bFieldVect(end)-_field.bFieldVect(beg)).mag();
// if the current interval is out of tolerance or too large, subdivide it    
    if ( test   >= _divTolerance || 
	 deltas >= divStepCeiling() ){
      //	 change >= _divTolerance ){
      divideRange(traj, slo, smid, posList);
      divideRange(traj, smid, sup, posList);
      return;
    }
  }
// otherwise, append these to the list and end
// append is OK to use here instead of insert since the
// above recursion guarantees the points will be entered in
// monotonic order.  If the above algorithm is changed, please
// make sure that this assumption still holds, or convert to
// using insert().
  posList.push_back(slo);
  posList.push_back(sup);
  return;
}









