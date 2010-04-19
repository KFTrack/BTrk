// ----------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPieceLineTraj.cc 524 2010-01-15 12:09:29Z stroili $
//
//  Description:
//  Continuous piecewise linear trajectory class
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 2/26/97
//-----------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrajGeom/TrkPieceLineTraj.hh"
#include "CLHEP/Geometry/HepPoint.h"
#include "BbrGeom/TrkGeomTrajVisitor.hh"

TrkPieceLineTraj::TrkPieceLineTraj(const TrkLineTraj& linetraj) :
  TrkPieceTraj(linetraj)
{
}

TrkPieceLineTraj::TrkPieceLineTraj(const HepPoint& point1, 
				   const HepPoint& point2) : 
  TrkPieceTraj(TrkLineTraj(point1,point2))
{
}

TrkPieceLineTraj::TrkPieceLineTraj(const TrkPieceLineTraj& other) :
  TrkPieceTraj(other)
{
}

TrkPieceLineTraj::~TrkPieceLineTraj()
{
}

TrkPieceLineTraj*
TrkPieceLineTraj::clone() const
{
  return new TrkPieceLineTraj(*this);
}

//  The full intellectual content of this class is the next operator
void
TrkPieceLineTraj::append(const HepPoint& hpoint)
{
//
//  Get the endpoint of the last trajectory.
//
  Trajectory* lasttraj = _traj.back().second;
  HepPoint lastpoint = lasttraj->position(lasttraj->hiRange());
//
//  Draw the line between this and the next point as the new line traj
//
  TrkLineTraj newtraj = TrkLineTraj(lastpoint,hpoint);
  TrkPieceTraj::append(newtraj);
}

void
TrkPieceLineTraj::accept(TrkGeomTrajVisitor& visitor) const
{
  visitor.visitPieceLine(this);
}
