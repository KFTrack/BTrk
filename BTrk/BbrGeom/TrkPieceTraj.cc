// ----------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPieceTraj.cc 491 2010-01-13 16:59:16Z stroili $
//
//  Description:
//  piecewise trajectory base class
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 2/26/97
//-----------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/TrkPieceTraj.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/BaBar/ErrLog.hh"
using std::cout;
using std::endl;
using std::ends;
using std::ostream;
using namespace CLHEP;
static const double STEPEPSILON = 1.0e-6; // .1 micron step

struct ordering {
        bool operator()(const std::pair<double,TrkGeomTraj*>& x,
                        const std::pair<double,TrkGeomTraj*>& y) const
        { return x.first < y.first; }
};

//- default constructor.  Make the vectors short
TrkPieceTraj::TrkPieceTraj() :
        TrkGeomTraj(-99999., 99999.)
{
  _traj.reserve(4);
}

//- construct from an existing trajectory
TrkPieceTraj::TrkPieceTraj(const TrkGeomTraj& seed) :
  TrkGeomTraj(0,seed.range())
{
  _traj.reserve(4);
  // intial global->local mapping has no offset
  _traj.push_back(std::make_pair(double(0), seed.clone() ) );
}

//- construct from an existing TrkPieceTraj
TrkPieceTraj::TrkPieceTraj(const TrkPieceTraj& other) :
  TrkGeomTraj(other.lowRange(),other.hiRange())
{
// deep-copy all the trajectory pieces
   _traj.reserve(other._traj.size());
   typedef TrkGeomTrajPV::const_iterator iter_t;
   for(iter_t i=other._traj.begin();i!=other._traj.end();++i) {
     _traj.push_back(std::make_pair( i->first, i->second->clone() ));
   }
}

//- assignment operator
TrkPieceTraj&
TrkPieceTraj::operator =(const TrkPieceTraj& other)
{
  flightrange[0] = other.flightrange[0];
  flightrange[1] = other.flightrange[1];
  deleteTraj();
//
// deep-copy all the trajectory pieces
//
  typedef TrkGeomTrajPV::const_iterator iter_t;
  _traj.reserve(other._traj.size());
  for(iter_t i=other._traj.begin();i!=other._traj.end();++i) {
    _traj.push_back(std::make_pair( i->first, i->second->clone() ));
  }
  return *this;
}

TrkPieceTraj::~TrkPieceTraj()
{
   deleteTraj();
}

void
TrkPieceTraj::deleteTraj()
{
  for (TrkGeomTrajPV::iterator i = _traj.begin();i!=_traj.end();i++) {
    delete i->second;
  }
  for( std::vector<TrkGeomTraj*>::iterator itrj = _removed.begin();itrj!= _removed.end();itrj++){
    delete *itrj;
  }
}

// append a trajectory
void
TrkPieceTraj::append(const TrkGeomTraj& nexttraj)
{
  _traj.push_back( std::make_pair( flightrange[1], nexttraj.clone() ));
  flightrange[1] += nexttraj.range();
}

// trajectory index from flight range
TrkPieceTraj::TrajIter
TrkPieceTraj::trajIndex(double flightdist,double& localflight) const
{
//  Return the appropriate end if the requested flightdistance is outside the
//  global range
  TrajIter i = flightdist < lowRange()  ? _traj.begin() :
               flightdist > hiRange()  ? _traj.end()-1 :
                 std::upper_bound(_traj.begin(),
                                  _traj.end(),
                                  std::make_pair(flightdist,(TrkGeomTraj*)0),
                                  ordering())-1;
#if 0
  if ( i->first > flightdist || flightdist > i->first+i->second->range()) {
   cout << "requested flightdist " << flightdist << " returning iterator to interval "
        << i->first << " -> " << i->first+i->second->range() << endl;
   printAll(cout);
  }
#endif
  localflight = localDist(i,flightdist);
  return i;
}

HepPoint
TrkPieceTraj::position(double flightdist) const
{
//
//  First, find the right trajectory piece, then give the local position
//
  double localflight(0.0);
  return trajIndex(flightdist,localflight)->second->position(localflight);
}

Hep3Vector
TrkPieceTraj::direction(double flightdist) const
{
//
//  First, find the right trajectory piece, then give the local direction
//
  double localflight(0.0);
  return trajIndex(flightdist,localflight)->second->direction(localflight);
}

double
TrkPieceTraj::curvature(double flightdist) const
{
//
//  First, find the right trajectory piece, then give the local curvature
//
  double localflight(0.0);
  return trajIndex(flightdist,localflight)->second->curvature(localflight);
}

Hep3Vector
TrkPieceTraj::delDirect(double flightdist) const
{
//
//  First, find the right trajectory piece, then give the local value
//
  double localflight(0.0);
  return trajIndex(flightdist,localflight)->second->delDirect(localflight);
}

void
TrkPieceTraj::getInfo(double flightdist,HepPoint& point,Hep3Vector& dir) const
{
//
//  First, find the right trajectory piece, then call the local function
//
  double localflight(0.0);
  trajIndex(flightdist,localflight)->second->getInfo(localflight,point,dir);
}

void
TrkPieceTraj::getInfo(double flightdist,HepPoint& point,Hep3Vector& dir,
                                        Hep3Vector& deldirect) const
{
//
//  First, find the right trajectory piece, then call the local function
//
  double localflight(0.0);
  trajIndex(flightdist,localflight)->second->getInfo(localflight,
                                                     point,dir,deldirect);
}

double
TrkPieceTraj::distTo1stError(double flightdist,double tol, int dir) const
{
//
//  First, find the right trajectory piece
//
  double localflight(0.0);
  TrajIter i = trajIndex(flightdist,localflight);
//
//  Ask the local piece for it's dist, and take the minimum of this or the
//  distance to the next trajectory piece. Take care of boundaries.
//
  double localdist = i->second->distTo1stError(localflight,tol,dir);
  double dist = localdist;
  if (dir > 0) {
    if( i+1 != _traj.end()) {
            double endofcurrent = i->first+i->second->range();
            dist = std::min(dist,endofcurrent-flightdist) + STEPEPSILON;
    }
  } else {
    if( i != _traj.begin()) {
            double beginofcurrent = i->first;
            dist = std::min(dist,flightdist - beginofcurrent) + STEPEPSILON;
    }
  }
  return dist;
}

double
TrkPieceTraj::distTo2ndError(double flightdist,double tol,int dir) const
{
//
//  First, find the right trajectory piece
//
  double localflight(0.0);
  TrajIter i = trajIndex(flightdist,localflight);
//
//  Ask the local piece for it's dist, and take the minimum of this or the
//  distance to the next trajectory piece. Take care of boundaries.
//
  double localdist = i->second->distTo2ndError(localflight,tol,dir);
  double dist = localdist;
  if (dir > 0) {
    if( i+1 != _traj.end()) {
            double endofcurrent = i->first+i->second->range();
            dist = std::min(dist,endofcurrent-flightdist) + STEPEPSILON;
    }
  } else {
    if( i != _traj.begin()) {
            double beginofcurrent = i->first;
            dist = std::min(dist,flightdist - beginofcurrent) + STEPEPSILON;
    }
  }
  return dist;
}


const TrkGeomTraj*
TrkPieceTraj::localTrajectory(double flightdist) const
{
//  Find and return the right trajectory piece
  double localflight(0);
  return trajIndex(flightdist,localflight)->second;
}
//
//  Print functions
//
void
TrkPieceTraj::print(ostream& os) const
{
  os << "TrkPieceTraj has " << _traj.size() << " pieces "
    << ", total flight range of " << hiRange();
}

void
TrkPieceTraj::printAll(ostream& os) const
{
  os << "TrkPieceTraj has " << _traj.size() << " pieces "
    << ", total flight range of " << hiRange() << endl;
  for(TrajIter i=_traj.begin();i!=_traj.end();++i){
    TrkGeomTraj* t = i->second;
    double lo = t->lowRange();
    double hi = t->hiRange();
    os << "Piece " << i-_traj.begin() << " has global range from "
       << i->first << " to " << i->first+t->range()
       << " and local range from " <<  lo << " to " << hi << endl;
    HepPoint start = t->position(lo);
    HepPoint end = t->position(hi);
    os << "Piece " << i-_traj.begin() << " starts at point "
       << start.x() <<","<< start.y()<<","<< start.z() 
       << " and ends at point "
       << end.x()<<"," << end.y()<<"," << end.z() 
       << endl;
  }
}
//
//  Re-size the last trajectory.  This can't extend into the
//  previous trajectory.
//

bool
TrkPieceTraj::resize(double locallen)
{
  TrajIter i = _traj.end()-1;
  double globallen = globalDist(i,locallen);
  bool goodlen = ( globallen > i->first );
  if(goodlen){
    flightrange[1] = globallen;
    double newrange[2];
    newrange[0] = i->second->lowRange();
    newrange[1] = locallen;
    i->second->setFlightRange(newrange);
  }
  return goodlen;
}
//
//  Range functions
//
void
TrkPieceTraj::setFlightRange(double newrange[2])
{
  double locdist;
  double lrange[2];
//
//  Check for pathological cases
//
  if( newrange[1] > newrange[0] &&
      newrange[0] < flightrange[1] &&
      newrange[1] > flightrange[0] ) {
//
//  Reset the lower range.  Delete trajectory pieces till we're in range.
//  This will leave the trajectory with a starting range different from
//  0, but it'll still be self-consistent
//
    TrajIter pivot = trajIndex(newrange[0],locdist);
    for (TrajIter i=_traj.begin();i!=pivot;++i) {
      _removed.push_back(i->second);
    }
    _traj.erase(_traj.begin(),_traj.begin()+(pivot-_traj.begin()));
//  reset the global range
    flightrange[0] = newrange[0];
    _traj.front().first = flightrange[0];
//  reset the local traj piece range
    lrange[0] = locdist;
    lrange[1] = _traj.front().second->hiRange();
    _traj.front().second->setFlightRange(lrange);
//
//  Same for the upper range
//
    pivot = trajIndex(newrange[1],locdist);
    for (TrajIter i = pivot+1; i!= _traj.end();++i) {
      _removed.push_back(i->second);
    }
    _traj.erase(_traj.begin()+(pivot+1-_traj.begin()),_traj.end());
    flightrange[1] = newrange[1];
    lrange[0] = _traj.back().second->lowRange();
    lrange[1] = locdist;
    _traj.back().second->setFlightRange(lrange);
  } else if(newrange[1] < newrange[0] ) {
    ErrMsg(error) << "TrkPieceTraj: cannot resize -- invalid range" << endmsg;
  } else {
//
//  Here the new range is completely discontinuous from the
//  old.  We'll just take the first (or last) traj piece,
//  clear out everything else, and reset the ranges.
//
    if( newrange[0] > flightrange[1]){
      for (TrajIter i = _traj.begin();i!=_traj.end()-1;++i) _removed.push_back(i->second);
      _traj.erase(_traj.begin(),_traj.end()-1);
    } else {
      for (TrajIter i = _traj.begin()+1; i!= _traj.end();++i) _removed.push_back(i->second);
      _traj.erase(_traj.begin()+1,_traj.end());
    }
    lrange[0] = localDist(_traj.begin(),newrange[0]);
    lrange[1] = localDist(_traj.begin(),newrange[1]);
    flightrange[0] = newrange[0];
    flightrange[1] = newrange[1];
    _traj.front().second->setFlightRange(lrange);
  }
}

const std::vector<std::pair<double,TrkGeomTraj*> >&
TrkPieceTraj::trajList() const
{
  return _traj;
}
