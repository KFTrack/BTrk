// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalEndSite.cc,v 1.20 2004/05/03 22:25:21 brownd Exp $
//
//  Description:
//  Trivial implementation of a KalSite to serve as the end of a processing
//  chain.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 7/2/97
//------------------------------------------------------------------------------
//
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalEndSite.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/TrkBase/TrkParams.hh"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
//
//  end site constructor
//
KalEndSite::KalEndSite(const TrkDifPieceTraj* reftraj,double len,trkDirection tdir,
		       double smearfactor,bool diagonly) :
  KalSite(endSite)
{
//  Set the trajectory
  setTraj(reftraj,len);
// construct the parameters
  params(tdir) = KalParams(*(localTrajectory()->parameters()));
// diagonalize if requested
  if(diagonly)params(tdir).diagonalize();
// smear the covariance
  params(tdir) *= smearfactor;
  setFit(tdir);
}

KalEndSite::KalEndSite(const KalParams& param,
		       const TrkDifPieceTraj* reftraj,
		       double len, trkDirection tdir,
		       double smearfactor,bool diagonly) :
  KalSite(endSite)
{
//  Set the trajectory
  setTraj(reftraj,len);
// construct the parameters
  params(tdir) = param;
// diagonalize if requested
  if(diagonly)params(tdir).diagonalize();
// smear the covariance
  params(tdir) *= smearfactor;
  setFit(tdir);
}

KalEndSite::KalEndSite(const KalEndSite& other) :
  KalSite(other)
{;}

KalSite*
KalEndSite::clone(const KalRep* ) const {
  return new KalEndSite(*this);
}

KalEndSite::KalEndSite() : KalSite(endSite)
{}

bool
KalEndSite::update(const TrkDifPieceTraj* reftraj,double len) {
  setTraj(reftraj,len);
  return true;
}
