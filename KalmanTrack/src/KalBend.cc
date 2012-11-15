// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalBend.cc,v 1.28 2006/04/24 18:53:05 brownd Exp $
//
//  Description: KalBend
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/30/98
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalBend.hh"
#include "TrkBase/TrkDifPieceTraj.hh"
#include "TrkBase/TrkParams.hh"
using std::endl;
using std::ostream;

KalBend::KalBend(const BFieldIntegrator& integrator,
		 const TrkDifPieceTraj* reftraj,BFieldIntRange const& range,
		 double momentum, int charge) :
  KalSite(KalSite::bendSite),_integrator(integrator),_range(range),_momentum(momentum),
  _charge(charge)
{
// update the cache
  updateCache(reftraj);
}
//
KalBend::KalBend(const KalBend& other) :
  KalSite(other),
  _integrator(other._integrator),
  _range(other._range),
  _momentum(other._momentum),
  _charge(other._charge),
  _delmom(other._delmom),
  _transport(other._transport),
  _momhat(other._momhat),
  _thetahat(other._thetahat),
  _phihat(other._phihat)
{
}
//
KalBend::~KalBend() {;}
//
KalBend*
KalBend::clone(const KalRep* newrep) const {
  KalBend* newbend = new KalBend(*this);
  newbend->setTraj(newrep->referenceTraj(),globalLength());
  return newbend;
}
//
// process; this is identical to KalMaterial
bool
KalBend::process(const KalSite* prevsite,trkDirection idir) {
  return processParams(prevsite,idir,_transport);
}
// update for a new trajectory
bool
KalBend::update(const TrkDifPieceTraj* newtraj,double newmom) {
  _momentum = newmom;
  updateCache(newtraj);
  return true;
}
//
void
KalBend::updateCache(const TrkDifPieceTraj* reftraj) {
// Set the trajectory parameters
// only update if the parameters have changed
  if(setTraj(reftraj,midpoint())){
// define the direction vectors
    _momhat = localTrajectory()->direction(localLength());
    _phihat = Hep3Vector(-_momhat.y(),_momhat.x(),0.0).unit();
//  _thetahat = _momhat.cross(_phihat);
    _thetahat = _momhat.cross(_phihat).unit();
// integrate the bfield over the range
    _delmom = _integrator.deltaMomentum(reftraj,_range);
// find the dependence of the parameters on deflection in the 2 views
    HepMatrix t1deflect = localTrajectory()->derivDeflect(localLength(),theta1);
    HepMatrix t2deflect = localTrajectory()->derivDeflect(localLength(),theta2);
// define the transport vector by applying these.  Note that no smearing
// is added.  The delta* functions take account of the charge-dependence
// of the bend corrections.
    double dtheta = deltaTheta();
    double dphi = deltaPhi();
    HepVector tvect(t1deflect*dtheta + t2deflect*dphi);
    _transport = KalParams(-tvect,HepSymMatrix(tvect.num_row(),0));
  }
// reset the site
  reset();
}

void
KalBend::printAll(ostream& os) const {
// call down to the base class first
  os << "Bend ";
  KalSite::printAll(os);
// print specific information.
  os << "Bend integral from " << _range._slo << " to "
     << _range._shi << endl;
  os << "Momentum change over integral = " << _delmom << endl;
  os << "Parameter transport over integral = " << _transport.parameterVector() << endl;
}

void
KalBend::invert() {
// Invert the range
  _range.invert();
// charge is inverted too!
  _charge =-_charge;
// call down to base
  KalSite::invert();
}
