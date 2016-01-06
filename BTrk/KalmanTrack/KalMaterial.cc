// File and Version Information:
//      $Id: KalMaterial.cc,v 1.70 2008/05/22 23:44:46 brownd Exp $
//
//  Description:
//  Class to describe a kalman filter material interaction site.
//  This models all material effects which do not change the particle
//  count (essentially multiple scattering and energy loss).
//
// Copyright Information:
//      Copyright (C) 1996      Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 12/18/96
//------------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include <math.h>
#include "BTrk/KalmanTrack/KalMaterial.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/DetectorModel/DetElem.hh"
#include "BTrk/DetectorModel/DetMaterial.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/TrkBase/TrkParams.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include <assert.h>
using std::endl;
using std::ostream;
using namespace CLHEP;
//

KalMaterial::KalMaterial(const DetIntersection& dinter,const TrkDifPieceTraj* reftraj,
    double momentum,TrkParticle const& tpart):
  KalSite(matSite),
  _tpart(tpart),
  _dinter(dinter),
  _momentum(momentum),
  _active(true),
  _pfract(0),_pfractrms(0),_deflectrms(0)
{
  assert(_dinter.trajet == (Trajectory*)reftraj); // check consistency
  updateCache(reftraj); // set the local cache
}
//
//  copy constructor
//
KalMaterial::KalMaterial(const KalMaterial& other) :
  KalSite(other),
  _tpart(other._tpart),
  _dinter(other._dinter),
  _transport(other._transport),
  _momentum(other._momentum),
  _active(other._active),
  _scatter(other._scatter),
  _eloss(other._eloss),
  _pderiv(other._pderiv),
  _pfract(other._pfract),
  _pfractrms(other._pfractrms),
  _deflectrms(other._deflectrms)
{}
// override setTraj function
bool
KalMaterial::setTraj(const TrkDifPieceTraj* reftraj,double globlen) {
  // Move detIntersection to new trajectory
  _dinter.trajet = reftraj;
  _dinter.pathlen = globlen;
  // call down to the base class function
  return KalSite::setTraj(reftraj,globlen);
}
//
KalMaterial::~KalMaterial(){;}
//
//  Update the site for a new intersection.
//
bool
KalMaterial::update(const TrkDifPieceTraj* newtraj,double newmom) {
  // update intersection
  _active = _dinter.delem->reIntersect(newtraj,_dinter);
  _momentum = newmom;
  updateCache(newtraj);
  return true;
}
//
//  process the effect of this materialing on fit result
//
bool
KalMaterial::process(const KalSite* prevsite,trkDirection idir){
  // if the change in momentum is small, use the linear approximation
  static const double largedp(0.01); // 1% cutoff
  bool status(false);
  if(_pfract < largedp)
    status =  processParams(prevsite,idir,_transport);
  else {
    // do the full non-linear propagation of momentum change
    invalidateSite(idir);
    // previous site must have been processed
    if(prevsite != 0 && prevsite->hasFit(idir)){
      if(isActive()){
        // copy over the previous site's parameters 
        copyParams(prevsite,idir);
        if(params(idir).matrixOK()){
          // process parameters analytically for this change in momentum.  Account for the convention
          // that the site's momentum is outside the site
          if(idir == trkIn)
            processDeltaP(params(idir),_pfract);
          else
            processDeltaP(params(idir),-_pfract/(1+_pfract));
          // add in the (linear) error
          params(idir).covarianceMatrix() += _transport.covarianceMatrix();
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
  }
  return status;
}
//
//  access
//
void
KalMaterial::printAll(ostream& os) const {
  os << "Material ";
  KalSite::printAll(os);
  os << "Transport vector, covariance = " << endl;
  _transport.printAll(os);
  os << "Intersection with "; detElem()->print(os);
  os << " Intersection flightlength, length = " << _dinter.pathlen << " "
    << _dinter.pathLength() << endl;
  os << "Fractional momentum loss =" << _pfract << " , RMS = " << _pfractrms << endl;
}
//
//  Update the cache when the underlying trajectory changes
//
void
KalMaterial::updateCache(const TrkDifPieceTraj* reftraj){
  // Set the trajectory parameters
  // only update derivatives if the parameters have changed
  if(setTraj(reftraj,_dinter.pathlen)){
    //  Get the material interaction information from the element.  By convention the momentum
    // provided is outwards of this site
    detElem()->materialInfo(_dinter,_momentum,_tpart,_deflectrms,_pfractrms,_pfract,trkIn);
    //  compute the (orthogonal) effects of scattering in the 2 directions
    //  on the covariance matrix.
    HepMatrix t1deflect = localTrajectory()->derivDeflect(localLength(),theta1);
    HepMatrix t2deflect = localTrajectory()->derivDeflect(localLength(),theta2);
    // energy loss effects are uncorrelated with scattering to good approximation
    _pderiv = localTrajectory()->derivPFract(localLength());
    // Matrices to propagate effects into parameter space
    static const HepSymMatrix unit(1,1);
    HepSymMatrix eunit = unit*pow(_pfractrms,2);
    HepSymMatrix sunit = unit*pow(_deflectrms,2);
    _eloss = eunit.similarity(_pderiv);
    // if this is a thick material, compute the transverse displacement effect as well
    static const double thickdeflect(0.01); // 100 microns average scattering in the material
    if(_deflectrms*_dinter.pathLength() < thickdeflect) {
      _scatter = sunit.similarity(t1deflect) + sunit.similarity(t2deflect);
    } else {
      // displacement effects: these are highly correlated with directional scattering
      HepMatrix t1displace = localTrajectory()->derivDisplace(localLength(),theta1);
      HepMatrix t2displace = localTrajectory()->derivDisplace(localLength(),theta2);
      // add these together according to formulas in PDG section 23.3, scaled by the path
      // length (converts angle to transverse displacement)
      // add scattering effects: first correlated angle and displacement (2 views are uncorrelated)
      static const double invsqrt3(1.0/sqrt(3.0));
      double cfac = 0.5*_dinter.pathLength();
      double ufac = cfac*invsqrt3;
      HepMatrix t1cor = t1deflect + cfac*t1displace;
      HepMatrix t2cor = t2deflect + cfac*t2displace;
      // then uncorrelated displacement effects
      HepMatrix t1uncor = ufac*t1displace;
      HepMatrix t2uncor = ufac*t2displace;
      // now create the effect matrices by expanding these into the parameter space
      _scatter = sunit.similarity(t1cor) + sunit.similarity(t2cor) + 
        sunit.similarity(t1uncor) + sunit.similarity(t2uncor);
    }
  }
  // compute the linear parameter transport
  _transport = KalParams(_pderiv*_pfract,_scatter+_eloss);
  // reset the site flags
  reset();
}

double
KalMaterial::radiationFraction() const {
  // get the material from the intersection
  return detElem()->material(_dinter).radiationFraction(_dinter.pathLength());
}

void
KalMaterial::setActivity(bool active) {
  if(active != _active){
    _active = active;
    setFit(trkIn, false);
    setFit(trkOut, false);
  }
}

double
KalMaterial::energyChange(trkDirection tdir) const {
  double energy = DetMaterial::particleEnergy(_momentum,_tpart);
  return momentumChange(tdir)*_momentum/energy;
}

double
KalMaterial::momentumChange(trkDirection tdir) const {
  return tdir == trkIn ? _pfract*_momentum : -_pfract*_momentum;
}

void
KalMaterial::invert() {
  // invert intersection fltlen
  _dinter.pathlen = -_dinter.pathlen;
  double temp = -_dinter.pathrange[0];
  _dinter.pathrange[0] = -_dinter.pathrange[1] ;
  _dinter.pathrange[1] = temp;
  // call down to base class
  KalSite::invert();
}
