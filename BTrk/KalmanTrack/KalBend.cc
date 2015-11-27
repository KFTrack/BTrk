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
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BField/BField.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalBend.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/TrkBase/TrkParams.hh"
using std::endl;
using std::ostream;
using namespace CLHEP;

double KalBend::_efac(0.0);
double KalBend::_berr(0.0);

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
KalBend::~KalBend() {;}
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
// deflection
    HepMatrix t1deflect = localTrajectory()->derivDeflect(localLength(),theta1);
    HepMatrix t2deflect = localTrajectory()->derivDeflect(localLength(),theta2);
// displacement
    HepMatrix t1displace = localTrajectory()->derivDisplace(localLength(),theta1);
    HepMatrix t2displace = localTrajectory()->derivDisplace(localLength(),theta2);
 // define the transport vector by applying these.  Note that no smearing
// is added.  The delta* functions take account of the charge-dependence
// of the bend corrections.
// angular deflection over this range
    double dtheta = deltaTheta();
    double dphi = deltaPhi();
    HepVector tvect(t1deflect*dtheta + t2deflect*dphi);
// displacement over this range
    HepVector dvect(t1displace*dtheta*_range.range() + t2displace*dphi*_range.range());
// compute the dimensionless error matrix.  The first term is from the uncertainty in the
// field measurement and interpolation, the second from the tracjectory
// notice I use the total momentum, not the transverse, as the field errors can be in any direction
    double berr2 = pow(BField::mmTeslaToMeVc*_range.range()*_berr/_momentum,2);
// the correction error scales as the sqrt of the correction, so that it is additive
    double eerr2 = _efac*_delmom.mag()/_momentum;
    static const HepSymMatrix unit(1,1);
    HepSymMatrix emat = unit*(berr2 + eerr2);
    static const double invsqrt3(1.0/sqrt(3.0));
// assume a 50% correlation between the deflection and displacement errors.  This models a random walk.    
    const double cfac = 0.5*_range.range();
    const double ufac = cfac*invsqrt3;
    HepMatrix t1cor = t1deflect + cfac*t1displace;
    HepMatrix t2cor = t2deflect + cfac*t2displace;
    HepMatrix t1uncor = ufac*t1displace;
    HepMatrix t2uncor = ufac*t2displace;
// correlated terms add linearly (amplitude), uncorrelated add in quadrature
    HepSymMatrix scatter = emat.similarity(t1cor) + emat.similarity(t2cor) + emat.similarity(t1uncor) + emat.similarity(t2uncor);
    _transport = KalParams(-tvect-dvect,scatter);
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
