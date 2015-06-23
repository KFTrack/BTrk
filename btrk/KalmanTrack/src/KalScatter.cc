// File and Version Information:
//      $Id: KalScatter.cc,v 1.3 2006/04/24 18:53:07 brownd Exp $
//
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 12/18/96
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include <math.h>
#include "KalmanTrack/KalScatter.hh"
#include "KalmanTrack/KalMaterial.hh"
#include "KalmanTrack/KalRep.hh"
#include "DetectorModel/DetElem.hh"
#include "DetectorModel/DetMaterial.hh"
#include "TrkBase/TrkDifPieceTraj.hh"
#include "TrkBase/TrkParams.hh"
#include <assert.h>
using std::endl;
using std::ostream;

KalScatter::KalScatter(const KalMaterial& km,double factor) :
  KalSite(scatterSite),
  _deflectrms(km.deflectRMS()*factor),
  _pfractrms(km.momFractionRMS()*factor)
{
  updateCache(static_cast<const TrkDifPieceTraj*>(km.detIntersection().trajet),km.globalLength());
}

//
//  copy constructor
//
KalScatter::KalScatter(const KalScatter& other) :
  KalSite(other),
  _transport(other._transport),
  _deflectrms(other._deflectrms),
  _pfractrms(other._pfractrms),
  _scatter(other._scatter)
{}
// clone function
KalScatter*
KalScatter::clone(const KalRep* krep) const {
  KalScatter* newscat = new KalScatter(*this);
  return newscat;
}
//
KalScatter::~KalScatter(){;}
//
//  Update the site for a new intersection.
//
bool
KalScatter::update(const TrkDifPieceTraj* newtraj,double newmom) {
  updateCache(newtraj,globalLength());
  return true;
}
//
//  process the effect of this scattering on fit result
//
bool
KalScatter::process(const KalSite* prevsite,trkDirection idir){
  return processParams(prevsite,idir,_transport);
}
//
//  access
//
void
KalScatter::printAll(ostream& os) const {
  os << "Scatter ";
  KalSite::printAll(os);
  os << "Transport vector, covariance = " << endl;
  _transport.printAll(os);
}
//
//  Update the cache when the underlying trajectory changes
//
void
KalScatter::updateCache(const TrkDifPieceTraj* reftraj,double fltlen){
// Set the trajectory parameters
// only update derivatives if the parameters have changed
  if(setTraj(reftraj,fltlen)){
//  compute the (orthogonal) effects of scattering in the 2 directions
//  on the covariance matrix.  Only update these if the trajectory has changed
    HepMatrix t1deflect = localTrajectory()->derivDeflect(localLength(),theta1);
    HepMatrix t2deflect = localTrajectory()->derivDeflect(localLength(),theta2);
    HepMatrix pderiv = localTrajectory()->derivPFract(localLength());
// now create the effect matrices by expanding these into the parameter space
    HepSymMatrix unit(1,1);
    double drms = pow(_deflectrms,2);
    double prms = pow(_pfractrms,2);
    _scatter = unit.similarity(t1deflect)*drms + 
      unit.similarity(t2deflect)*drms + 
      unit.similarity(pderiv)*prms;
  }
  HepVector null(_scatter.num_row(),0); // null effect vector
  _transport = KalParams(null,_scatter);
// reset the site flags
  reset();
}
