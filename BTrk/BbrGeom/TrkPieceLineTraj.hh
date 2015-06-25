// ----------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPieceLineTraj.hh 524 2010-01-15 12:09:29Z stroili $
//
//  Description:
//  Continuous piecewise linear trajectory class
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 2/26/97
//-----------------------------------------------------------------------------
#ifndef TRKPIECELINETRAJ_HH
#define TRKPIECELINETRAJ_HH
#include "BTrk/BbrGeom/TrkPieceTraj.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
class HepPoint;

class TrkPieceLineTraj : public TrkPieceTraj
{
public:
//  construct from a line traj (or a pair of points)
  TrkPieceLineTraj( const TrkLineTraj&);
  TrkPieceLineTraj( const HepPoint&,const HepPoint&);
  TrkPieceLineTraj( const TrkPieceLineTraj&  );   // copy ctor
  ~TrkPieceLineTraj();

  TrkPieceLineTraj* clone() const;

  // Support Visitor pattern (see TrkGeomTraj.hh)
  void accept(TrkGeomTrajVisitor&) const;

//  add-a-point operator; this builds the line trajectories
  void append(const HepPoint&);
};
#endif
