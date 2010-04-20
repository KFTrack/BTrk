//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalContext.cc 103 2010-01-15 12:12:27Z stroili $
//
// Description:
//      class KalContext
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 3/15/97
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkEnv/KalContext.hh"
#include "AbsEnv/AbsEnv.hh"
#include "TrkEnv/TrkEnv.hh"

KalContext::KalContext() :
  _disttol(0.1), // spatial separation between trajs to continue iterating
  _maxiter(3), // maximum number of iteration steps
  _matsites(true), // use material intersections as parameter references
  _bends(true), // KalBend correction for non-homogeneous Bfield
  _smearfactor(1.0e8), // matrix smearing, should depend on P
  _sitethresh(0.2), // momentum fraction change to stop the track in a site
  _momthresh(0.5), // total momentum fraction change to stop the track
  _sitepfrac(0.01), // fractional momentum loss to trigger local parameter reference
  _sitedflct(0.01), // deflection to trigger local parameter reference
  _mingap(1.0e-4), // Minimum flight distance gap to create a new trajectory piece
  _trajbuff(0.001), // small buffer when appending trajectories
  _bintminstep(0.5), // BField integration parameters
  _bintmaxstep(5.0),
  _bintmaxfrac(0.1),
  _binttolerance(0.01),
  _bdivminstep(0.5), // BField track divider parameters
  _bdivmaxstep(5.0),
  _bdivmaxfrac(0.1),
  _bdivtolerance(0.01),
  _defpid(PdtPid::pion),
  _maxmomdiff(0.05),
  _stophots(false),
  _ambigflip(true),
  _momfac(0.0)
{
// max par diff in units chi^2 units: note trkOut=0, trkIn=1 !!!
  _maxpardif[0] = _maxpardif[1] = 1.0; // parameter pull difference for iteration convergence testing
// default volumes span the global tracking volume
  _volumes[trkIn] = TrkVolumeHandle::beampipe;
  _volumes[trkOut] = TrkVolumeHandle::dch;
// DOF requirements based on helix assumption
  _mindof[TrkEnums::xyView] = 3;
  _mindof[TrkEnums::zView] = 2;
  _mindof[TrkEnums::bothView] = 1;
}

KalContext::~KalContext(){}

const TrkVolume*
KalContext::trkVolume(trkDirection tdir) const {
// cast to the enum
  TrkVolumeHandle::trkvolumes vol = (TrkVolumeHandle::trkvolumes)_volumes[tdir];
  return gblEnv->getTrk()->findTrkVolume(vol);
}
    
