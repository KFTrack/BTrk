// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalEndSite.hh,v 1.21 2004/04/14 04:59:12 brownd Exp $
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
#ifndef KALENDSITE_HH
#define KALENDSITE_HH

#include "KalmanTrack/KalSite.hh"
#include "ErrLogger/ErrLog.hh"

class KalEndSite : public KalSite {
public:
//  Construct from piece traj
  KalEndSite(const TrkDifPieceTraj* ptraj,double fltlen,trkDirection,
	     double smearfactor=1.0,bool diagonly=false);
// Construct from parameters; piecetraj is just used for reference point
  KalEndSite(const KalParams&,
	     const TrkDifPieceTraj* ptraj,
	     double fltlen,trkDirection,
	     double smearfactor=1.0,bool diagonly=false);
// copy constructor
  KalEndSite(const KalEndSite&);
// clone operator
  KalSite* clone(const KalRep*) const;
  virtual ~KalEndSite() {;}
// allow a null constructor
  KalEndSite();
// Procesing
  bool process(const KalSite*,trkDirection)
    { ErrMsg(error) << "KalEndSite:Error: process function called " << endmsg;
    return false; }
  virtual bool update(const TrkDifPieceTraj*,double);
};

#endif
