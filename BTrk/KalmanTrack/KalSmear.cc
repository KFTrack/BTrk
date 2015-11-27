// File and Version Information:
//      $Id: KalSmear.cc,v 1.6 2006/04/24 18:53:08 brownd Exp $
//
//      class KalSmear
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 12/18/96
//------------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalSmear.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/TrkBase/TrkParams.hh"

#include "CLHEP/Matrix/SymMatrix.h"
#include <assert.h>
using std::endl;
using std::ostream;
using namespace CLHEP;

//
KalSmear::KalSmear(const HepSymMatrix& smear,
		   const TrkDifPieceTraj* ptraj,
		   double glen) :
  KalSite(smearSite)
{
// initialize the trajectory
  setTraj(ptraj,glen);
// setup the transport tensor
  HepVector pvec(smear.num_row(),0);
  _transport = KalParams(pvec,smear);
}
//
KalSmear::KalSmear(const TrkDifPieceTraj* ptraj,
		   double glen,
		   double smearfac) :
  KalSite(smearSite)
{
// initialize the trajectory
  setTraj(ptraj,glen);
// construct the parameters (0 transport vector)
  KalParams params(*(localTrajectory()->parameters()));
  HepVector pvec(params.covarianceMatrix().num_row(),0);
  _transport = KalParams(pvec,params.covarianceMatrix());
// smear the covariance
  _transport *= smearfac;
// 'diagonalize' it too
  _transport.diagonalize();
}
//
KalSmear::~KalSmear(){;}
//
//  Update the site for a new intersection.
//
bool
KalSmear::update(const TrkDifPieceTraj* newtraj,double) {
// nothing to do besides reset the base
  setTraj(newtraj,globalLength());
  reset();
  return true;
}
//
//  process the effect of this smearing on fit result
//
bool
KalSmear::process(const KalSite* prevsite,trkDirection idir){
  return processParams(prevsite,idir,_transport);
}
//
//  access
//
void
KalSmear::printAll(ostream& os) const {
  os << "KalSmear ";
  KalSite::printAll(os);
  os << "Transport vector, covariance = " << endl;
  _transport.printAll(os);
}
