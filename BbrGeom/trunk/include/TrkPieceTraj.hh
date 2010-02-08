// ----------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPieceTraj.hh 491 2010-01-13 16:59:16Z stroili $
//
//  Description:
//  Piecewise compound trajectory base class.  The useful constructors and the
//  'append' function are protected, insuring that they'll only be used by
//  subclasses.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 2/26/97
//-----------------------------------------------------------------------------
#ifndef TRKPIECETRAJ_HH
#define TRKPIECETRAJ_HH
#include "BbrGeom/TrkGeomTraj.hh"
#include <vector>
#include <utility>

#include <iosfwd>

class TrkPieceTraj : public TrkGeomTraj
{
public:
//
//  dummy constructor
//
  TrkPieceTraj();
// destructor; pure virtual, to make the class abstract
  virtual ~TrkPieceTraj();
// operators
  TrkPieceTraj& operator=(const TrkPieceTraj&);

// Normal trajectory functions: these are all fully implemented as call downs
// to the constituent trajectories.
  HepPoint position( double ) const;
  Hep3Vector direction( double ) const;
  virtual double curvature( double f = 0 ) const;
  Hep3Vector delDirect( double ) const;
  void  getInfo(double fltLen, HepPoint& , Hep3Vector& direction) const;
  void  getInfo(double fltLen, HepPoint& , Hep3Vector& direction,
                               Hep3Vector& delDirect) const;
  double distTo1stError(double s, double tol, int pathDir) const;
  double distTo2ndError(double s, double tol, int pathDir) const;

//  Piecewise-traj specific functions
  const TrkGeomTraj* localTrajectory(double) const;

//  Overwrite the range functions, to allow local ranges to be
//  updated as well
  void setFlightRange(double newrange[2]);

  void print(std::ostream& os) const;
  void printAll(std::ostream& os) const;

  // Provide the list of Trajs (to permit Visitors to function).
  const std::vector<std::pair<double,TrkGeomTraj*> >& trajList() const;

protected:
  //                            global fl., traj
  typedef std::vector<std::pair<double,TrkGeomTraj*> > TrkGeomTrajPV;
  typedef TrkGeomTrajPV::const_iterator TrajIter;
  TrkPieceTraj( const TrkGeomTraj&  ); // build from a seed
  TrkPieceTraj(const TrkPieceTraj&); // copy constructor
  void append(const TrkGeomTraj&); //  append-a-trajectory function
  bool resize(double locallen); // extend/shorten the last trajectory
  TrajIter trajIndex(double globalDist, double& localDist) const;
  inline double localDist(TrajIter index,double globdist) const;
  inline double globalDist(TrajIter index,double locdist) const;
// Data Members
  TrkGeomTrajPV _traj; // Ordered vector of trajs
  std::vector<TrkGeomTraj*> _removed; // vector of removed trajs.  We can't delete these until the whole object is deleted
  void deleteTraj();
};

double
TrkPieceTraj::localDist(TrajIter index,double globdist) const
{
  return index->second->lowRange() + ( globdist - index->first );
};

double
TrkPieceTraj::globalDist(TrajIter index,double locdist) const
{
  return index->first + ( locdist - index->second->lowRange() );
};
#endif
