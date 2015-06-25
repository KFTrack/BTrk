//--------------------------------------------------------------------------
//$Id: BFieldIntegrator.cc 830 2011-04-09 07:14:58Z brownd $
//Description:Class Implementation for |BFieldIntegrator|
//Author List:A. Snyder
//Copyright (C) 1998	SLAC
//------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "BTrk/BField/BField.hh"
#include "BTrk/BbrGeom/Trajectory.hh"
#include "BTrk/BField/BFieldIntegrator.hh"
#include <vector>
using namespace CLHEP;


BFieldIntegrator::BFieldIntegrator(const BField &f,BFieldIntConfig const& config) :
  _field(f),_config(config),_bNominal(0.0, 0.0, f.bFieldNominal())
{}


BFieldIntegrator::BFieldIntegrator(const BField &f,BFieldIntConfig const& config,const Hep3Vector& bnom) :
  _field(f),_config(config), _bNominal(bnom)
{}

BFieldIntegrator::~BFieldIntegrator()
{}

Hep3Vector
BFieldIntegrator::deltaMomentum(const Trajectory *traj, BFieldIntRange const& range) const
{
  static Hep3Vector null(0.0,0.0,0.0);
  Hep3Vector retval = null;
  HepPoint beg = traj->position(range._slo);
  HepPoint end = traj->position(range._shi);
  double deltas = range._shi - range._slo;
  if(deltas <  _config._intPathMin ) {
    HepPoint beg=traj->position(range._slo);
    HepPoint end=traj->position(range._shi);
    HepPoint mid=traj->position(range._smid);
    Hep3Vector bbeg(_field.bFieldVect(beg));
    Hep3Vector bend(_field.bFieldVect(end));
    Hep3Vector bmid(_field.bFieldVect(mid));
    double db = (bend -bbeg).mag();
    Hep3Vector bdiff = bmid - _bNominal;
// compute the change in momentum from this
    double dpmax = deltas * db * BField::mmTeslaToMeVc;
    if (dpmax < _config._intTolerance){
      Hep3Vector dir = traj->direction(range._smid) * deltas * BField::mmTeslaToMeVc;
      retval = dir.cross(bdiff);
    }
  } else if(deltas < _config._maxRange) {
    retval = deltaMomentum(traj, BFieldIntRange(range._slo, range._smid)) + deltaMomentum(traj, BFieldIntRange(range._smid, range._shi));
  } 
  return retval;
}

void
BFieldIntegrator::divideRange(const Trajectory *traj,BFieldIntRange const& range,
    std::vector<BFieldIntRange>& divrange) const
{
  bool divide(false);
  double deltas = range._shi - range._slo;
// trivial tests first
  if (deltas > _config._divPathMin && deltas < _config._maxRange) {
    divide = deltas > _config._divStepCeiling;
    if(!divide){
// calculate the fractional change in momentum over this range
      HepPoint beg=traj->position(range._slo);
      HepPoint end=traj->position(range._shi);
      HepPoint mid=traj->position(range._smid);
      Hep3Vector bbeg(_field.bFieldVect(beg));
      Hep3Vector bend(_field.bFieldVect(end));
      Hep3Vector bmid(_field.bFieldVect(mid));
      Hep3Vector bavg = (bbeg+bmid+bend)/3.0;
      // Find the largest of either the difference with the nominal, or the change over the range
      double db = (bend -bbeg).mag();
      double bdiff = (bavg - _bNominal).mag();
      double deltab = std::max(db,bdiff);
      // compute the change in momentum from this
      double dpmax = deltas * deltab * BField::mmTeslaToMeVc;
      divide = dpmax > _config._divTolerance;
    }
  }
  if(divide){
    divideRange(traj, BFieldIntRange(range._slo, range._smid), divrange);
    divideRange(traj, BFieldIntRange(range._smid, range._shi), divrange);
  } else {
    divrange.push_back(range);
  }
}

