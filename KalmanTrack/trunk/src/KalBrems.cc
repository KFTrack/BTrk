// File and Version Information:
//      $Id: KalBrems.cc,v 1.4 2006/04/27 17:59:14 brownd Exp $
//
//
// Copyright Information:
//	Copyright (C) 2006	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/7/2006
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "BaBar/Constants.hh"
#include <math.h>
#include "KalmanTrack/KalBrems.hh"
#include "KalmanTrack/KalRep.hh"
#include "AbsCalo/AbsRecoCalo.hh"
#include "CLHEP/Utilities/CLHEP.h"
#include "TrkBase/TrkDifPieceTraj.hh"
#include "TrkBase/TrkParams.hh"
#include <assert.h>
using std::endl;
using std::ostream;
using std::cout;
//


KalBrems::KalBrems(const TrkDifPieceTraj* reftraj,
		   const AbsRecoCalo* cluster,
		   double momentum,
		   PdtPid::PidType pid) :
  KalSite(bremsSite),
  _cluster(cluster),
  _partid(pid),
  _momentum(momentum),
  _active(pid == PdtPid::electron),
  _gamphi(0),
  _gamE(0.0),
  _gamEerr(0.0)
{
// set the local cache
  updateCache(reftraj);
}


//
//  copy constructor
//
KalBrems::KalBrems(const KalBrems& other) :
  KalSite(other),
  _cluster(other._cluster),
  _partid(other._partid),
  _momentum(other._momentum),
  _active(other._active),
  _gamphi(other._gamphi),
  _gamE(other._gamE),
  _gamEerr(other._gamEerr)
{}

// clone function
KalBrems*
KalBrems::clone(const KalRep* krep) const {
  KalBrems* newbrems = new KalBrems(*this);
// see if the PID has changed: if not, we can reset the traj
// and use the existing cache
  if(newbrems->_partid == krep->particleType()){
// set the trajectory to the new rep's reference trajectory
// For now, assume the length doesn't change
    newbrems->setTraj(krep->referenceTraj(),globalLength());
  } else {
//  If the PID has changed, we need to re-compute the cache
    newbrems->setPID(krep->particleType());
    newbrems->update(krep->referenceTraj(),krep->refMomentum());
  }
  return newbrems;
}
//
KalBrems::~KalBrems(){;}
//
//  Update the site for a new intersection.
//
bool
KalBrems::update(const TrkDifPieceTraj* newtraj,double newmom) {
  _momentum = newmom;
  updateCache(newtraj);
  return true;
}
//
//  process the effect of this materialing on fit result
//
bool
KalBrems::process(const KalSite* prevsite,trkDirection idir){
  bool status(false);
  invalidateSite(idir);
// previous site must have been processed
  if(prevsite != 0 && prevsite->hasFit(idir)){
    if(isActive()){
// copy over the previous site's parameters 
      copyParams(prevsite,idir);
      if(params(idir).matrixOK()){
// explicitly transport the parameters according to the 
        double pfract,pfractrms;
        if(idir == trkIn){
          pfract = clusterEnergy()/_momentum;
          pfractrms = clusterEnergyErr()/_momentum;
        } else{ 
          pfract = -clusterEnergy()/(clusterEnergy()+_momentum);
          pfractrms = clusterEnergyErr()/(clusterEnergy()+_momentum);
        }
        processDeltaP(params(idir),pfract);
// computer the error and add it in
        HepSymMatrix unit(1,1);
        unit *=sqr(pfractrms);
        HepSymMatrix gamerr = unit.similarity(_pderiv);
// update the sites errors
        params(idir).covarianceMatrix() += gamerr;
// set the site processed in this direction
        setFit(idir);
        status = true;
      }
    } else {
// copy the previous site
      copySite(prevsite,idir);
      if(validSite(idir)){
// set the site processed in this direction
        setFit(idir);
        status = true;
      }
    }
  }
  return status;
}

//
//  access
//
void
KalBrems::printAll(ostream& os) const {
  os << "Brems";
  KalSite::printAll(os);
  os << "Momentum = " << _momentum;
  os << " Cluster energy = " << clusterEnergy() << " +- " << clusterEnergyErr() << endl; 
}
//
//  Update the cache when the underlying trajectory changes
//
void
KalBrems::updateCache(const TrkDifPieceTraj* reftraj){
// find the tangent point
  double fltlen = tangent(reftraj);
// Set the trajectory parameters
  setTraj(reftraj,fltlen);
// energy loss effects derivatives
  _pderiv = localTrajectory()->derivPFract(localLength());
// compute the direction of the photon.  This is computed from the position on the
// trajectory to the cluster.  Eventually we should re-compute the tangent point in this
// function, instead of assuming the flight-length doesn't change. FIXME ******
  Hep3Vector gamdir = _cluster->direction(reftraj->position(globalLength()));
// extract phi
  _gamphi = gamdir.phi();
// set cluster energy
  _gamE = _cluster->energy(localTrajectory()->position(localLength()));
  _gamEerr = _cluster->energyErr(localTrajectory()->position(localLength()));
// reset the site flags
  reset();
}

void
KalBrems::setActivity(bool active) {
  if(active != _active){
    _active = active;
    setFit(trkIn, false);
    setFit(trkOut, false);
  }
}


void
KalBrems::setPID(PdtPid::PidType newtype){
  _partid = newtype;
  setActivity(newtype == PdtPid::electron);
}

double
KalBrems::clusterEnergy() const {
  return _gamE;
}

double
KalBrems::clusterEnergyErr() const {
  return _gamEerr;
}

double
KalBrems::momentumChange(trkDirection idir) const {
  if(idir == trkOut)
    return -clusterEnergy();
  else
    return clusterEnergy();
}

double
KalBrems::tangent(const TrkDifPieceTraj* reftraj) const {
// start with the origin 
  double tflt(0.0);
// setup iteration conditions
  static const double tol(0.1); //  tolerance on iteration
  static const unsigned maxiter(10); // max # of iterations
  unsigned iter(0);
  double delta(1.0e6);
// iterate till converged or we get tired
  bool diverge(false);
  while(fabs(delta) > tol && iter < maxiter &&  !diverge){
// compute the photon direction at this point
    HepPoint tpos = reftraj->position(tflt);
    Hep3Vector tdir = reftraj->direction(tflt);
    Hep3Vector ddir = reftraj->delDirect(tflt);
    Hep3Vector cdir = (_cluster->position(tpos)-tpos).unit();
// compute change in direction
    Hep3Vector deltadir = cdir - tdir;
    Hep3Vector ddirnorm = ddir.unit();
// compute change in flightlength this implies
    double newdelta = deltadir.dot(ddirnorm)/ddir.mag();
    diverge = fabs(newdelta/delta) > 2.0;
// update
    delta = newdelta;
    tflt += delta;
    iter++;
  }
  if(diverge) ErrMsg(error) << "Brems iterations diverging " << endmsg;
  if(iter == maxiter) ErrMsg(error) << "Brems iterations didn't converge" << endmsg;
  return tflt;
}
