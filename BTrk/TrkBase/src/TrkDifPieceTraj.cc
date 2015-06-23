// ---------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkDifPieceTraj.cc,v 1.52 2008/05/05 16:45:57 fransham Exp $
//
//  Description:
//  Differential form of the Piecewise compound trajectory base class.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/17/97
//----------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "BaBar/Constants.hh"
#include "TrkBase/TrkDifPieceTraj.hh"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/HepPoint.h"
#include "difAlgebra/DifPoint.hh"
#include "difAlgebra/DifVector.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkErrCode.hh"
#include "ErrLogger/ErrLog.hh"
#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include "BaBar/BbrCollectionUtils.hh"
using std::endl;
using std::ends;
using std::ostream;

static const double _STEPEPSILON = 1.0e-6; // .1 micron step
static const double _TOL = 1.0e-5; // poca tolerance
static const double _SMALLGAP = 0.1; // negligible trajectory gap

// trajectory comparison object
struct sortTrajs : std::binary_function<Trajectory*,Trajectory*,bool> {
  bool operator() ( Trajectory* const & a, Trajectory * const  & b){
    return a->lowRange() < b->lowRange();
  }
};

//
//  Size the vectors to a reasonable number of terms
//
TrkDifPieceTraj::TrkDifPieceTraj(const TrkSimpTraj& seed,const double lowlim,const double hilim) :
  TrkDifTraj(lowlim,hilim), _lastIndex(-1)
{
  // _localtraj.reserve(16); _globalrange.reserve(16);
  _localtraj.push_back((TrkSimpTraj*)seed.clone());
  assert(lowlim < hilim);
  _globalrange.push_back(lowlim);
  _globalrange.push_back(hilim);
// don't assume the local trajectory has the same range, but do respect it's starting point
  double locrange[2];
  locrange[0] = seed.lowRange();
  locrange[1] = hilim-lowlim+locrange[0];
  _localtraj.front()->setFlightRange(locrange);
}

TrkDifPieceTraj::TrkDifPieceTraj(TrkSimpTraj* seed,const double lowlim,const double hilim) :
  TrkDifTraj(lowlim,hilim), _lastIndex(-1)
{
  // _localtraj.reserve(16); _globalrange.reserve(16);
  assert(lowlim < hilim);
  _globalrange.push_back(lowlim);
  _globalrange.push_back(hilim);
// don't assume the local trajectory has the same range, but do respect it's starting point
  double locrange[2];
  locrange[0] = seed->lowRange();
  locrange[1] = hilim-lowlim+locrange[0];
  seed->setFlightRange(locrange);
  _localtraj.push_back(seed);
}

//
TrkDifPieceTraj::TrkDifPieceTraj(const TrkDifPieceTraj& other) :
  TrkDifTraj(other.lowRange(),other.hiRange()),
  _globalrange(other._globalrange),
  _lastIndex(other._lastIndex)
{
//
// deep-copy all the trajectory pieces
//
  // _localtraj.reserve(other._localtraj.size());
  typedef std::deque<TrkSimpTraj *>::const_iterator iter_t;
  iter_t end = other._localtraj.end();
  for(iter_t i=other._localtraj.begin();i!=end; ++i)
    _localtraj.push_back( (TrkSimpTraj*) (*i)->clone());
}

TrkDifPieceTraj::TrkDifPieceTraj(const std::vector<TrkSimpTraj*>& trajs)
{
  // _localtraj.reserve(trajs.size());_globalrange.reserve(trajs.size()+1);
// sort the trajectories.  This requires making a copy (ugh)
  std::vector<TrkSimpTraj*> strajs(trajs);
  std::sort(strajs.begin(),strajs.end(),sortTrajs());
// intialize
  TrkSimpTraj* prevtraj(0);
  for(std::vector<TrkSimpTraj*>::const_iterator itraj=strajs.begin();itraj!=strajs.end();itraj++){
    TrkSimpTraj* newtraj = (*itraj)->clone();
    assert(newtraj != 0);
    if(prevtraj != 0){
      TrkErrCode add;
      double gap;
      if(newtraj->lowRange() > prevtraj->lowRange())
        add = append(newtraj,gap);
      else
        add = prepend(newtraj,gap);
      if(add.failure()){
        ErrMsg(warning) << "construction from vector of trajs failed"
                        << add << endmsg;
        delete newtraj;
        break;
      }
    } else {
// Initial global flightlength is defined by the first trajectory
      _localtraj.push_back(newtraj);
      _globalrange.push_back(newtraj->lowRange());
      _globalrange.push_back(newtraj->hiRange());
// set the global range
      flightrange[0] = _globalrange.front();
      flightrange[1] = _globalrange.back();
    }
    prevtraj = _localtraj[itraj-strajs.begin()];
  }
}

TrkDifPieceTraj::TrkDifPieceTraj(std::vector<TrkSimpTraj*>& trajs,double gstart)
{
// set Initial global flightlength 
  _globalrange.push_back(gstart);
  for(std::vector<TrkSimpTraj*>::const_iterator itraj=trajs.begin();itraj!=trajs.end();itraj++){
    TrkSimpTraj* newtraj = *itraj;
    assert(newtraj != 0);
// check
    if(_localtraj.size()>0){
      HepPoint oldps = position(_globalrange.back());
      HepPoint newps = newtraj->position(newtraj->lowRange());
      if(oldps.distanceTo(newps) > 1e-3 ){
        ErrMsg(error) << "traj points don't match" << endmsg;
      }
    }
    _globalrange.push_back(_globalrange.back() + newtraj->hiRange() - newtraj->lowRange());
    _localtraj.push_back(newtraj);
  }
// set the global range
  flightrange[0] = _globalrange.front();
  flightrange[1] = _globalrange.back();
}

TrkDifPieceTraj::~TrkDifPieceTraj()
{
  std::for_each(_localtraj.begin(),
                _localtraj.end(),
                babar::Collection::DeleteObject());
}

void
TrkDifPieceTraj::getDFInfo(double flightdist, DifPoint& pos, DifVector& direction,
                           DifVector& delDirect) const
{
//
//  First, find the right trajectory piece, then give the local position
//
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(flightdist,localflight);
  loctraj->getDFInfo(localflight,pos,direction,delDirect);
}

void
TrkDifPieceTraj::getDFInfo2(double flightdist, DifPoint& pos,
                            DifVector& direction) const
{
//
//  First, find the right trajectory piece, then give the local position
//
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(flightdist,localflight);
  loctraj->getDFInfo2(localflight,pos,direction);
}
//- assignment operator
TrkDifPieceTraj&
TrkDifPieceTraj::operator =(const TrkDifPieceTraj& other)
{
  if (&other == this) return *this;
  flightrange[0] = other.flightrange[0];
  flightrange[1] = other.flightrange[1];
  _globalrange = other._globalrange;
//
// deep-copy all the trajectory pieces
//
  std::for_each(_localtraj.begin(),_localtraj.end(),babar::Collection::DeleteObject());
  _localtraj.clear(); // _localtraj.reserve(other._localtraj.size());
  typedef std::deque<TrkSimpTraj*>::const_iterator iter_t;
  for(iter_t i=other._localtraj.begin();i!=other._localtraj.end();++i)
    _localtraj.push_back((*i)->clone());
  return *this;
}
//  Append a new trajectory at a given global flight length.  This makes a linear
//  interpolation to set the range of the new trajectory to match the position
//  of the existing last piece.
//
const TrkErrCode&
TrkDifPieceTraj::append(double glen,TrkSimpTraj* nexttraj,double& gap)
{
  static TrkErrCode retval;
  retval = TrkErrCode(TrkErrCode::succeed);
// simple range checks
  if(glen >= lowRange()){
// resize the end if necessary
    int nremoved = resize(glen,trkOut);
    if(nremoved > 0) ErrMsg(warning) << "append removed "
                                     << nremoved << " old trajectories" << endmsg;
//  Compute POCA for this trajectory to the point at the end of the existing trajectory
    HepPoint endpoint = position(glen);
    HepPoint newend = nexttraj->position(nexttraj->lowRange());
    gap = newend.distanceTo(endpoint);
    double range[2];
    range[0] = nexttraj->lowRange();
    range[1] = nexttraj->hiRange();
    double delta = std::max(hiRange() - glen,range[1]-range[0]);
    if(gap > _TOL){ // Don't invoke POCA if the point is too close
      TrkPoca endpoca(*nexttraj,nexttraj->lowRange(),endpoint,_TOL);
      if(endpoca.status().failure()){
        retval = TrkErrCode(TrkErrCode::fail,7,"poca failure");
        return retval;
      }
      range[0] = endpoca.flt1();
      gap = endpoca.doca();
    }
// make sure the global range stays at least as long as before
    range[1] = range[0] + delta;
// resize it
    nexttraj->setFlightRange(range);
// append the trajectory
    _localtraj.push_back(nexttraj);
// append the new global range
    _globalrange.push_back(glen+delta);
    flightrange[1] = glen+delta;
  } else
    retval = TrkErrCode(TrkErrCode::fail,6," cannot append before minimum flight range");
  return retval;
}
// similair code for prepending
const TrkErrCode&
TrkDifPieceTraj::prepend(double glen,TrkSimpTraj* nexttraj,double& gap)
{
  static TrkErrCode retval;
  retval = TrkErrCode(TrkErrCode::succeed);
// simple range checks
  if(glen <= hiRange()){
// resize
    int nremoved = resize(glen,trkIn);
    if(nremoved > 0) ErrMsg(warning) <<  " prepend removed "
                                     << nremoved << " old trajectories" << endmsg;
//  Compute POCA for this trajectory to the point at the end of the existing trajectory
    HepPoint endpoint = position(glen);
    HepPoint newend = nexttraj->position(nexttraj->hiRange());
    gap = newend.distanceTo(endpoint);
    double range[2];
    range[0] = nexttraj->lowRange();
    range[1] = nexttraj->hiRange();
    double delta = std::max(glen -lowRange(),range[1]-range[0]);
    if(gap > _TOL){ // Don't invoke POCA if the point is too close
      TrkPoca endpoca(*nexttraj,nexttraj->hiRange(),endpoint,_TOL);
      if(endpoca.status().failure()){
        retval = TrkErrCode(TrkErrCode::fail,7,"poca failure");
        return retval;
      }
      range[1] = endpoca.flt1();
      gap = endpoca.doca();
    }
// make sure the global range stays at least as long as before
    range[0] = range[1] - delta;
// resize it
    nexttraj->setFlightRange(range);
// prepend the trajectory
    _localtraj.push_front(nexttraj);
// prepend the new global range
    _globalrange.push_front(glen-delta);
    flightrange[0] = glen-delta;
  } else
    retval = TrkErrCode(TrkErrCode::fail,6," cannot prepend after maximum flight range");
  return retval;
}

const TrkErrCode&
TrkDifPieceTraj::append(double glen,const TrkSimpTraj& nexttraj,double& gap)
{
  return append(glen,nexttraj.clone(),gap);
}

const TrkErrCode&
TrkDifPieceTraj::prepend(double glen,const TrkSimpTraj& nexttraj,double& gap)
{
  return prepend(glen,nexttraj.clone(),gap);
}

// append an entire TrkDifTraj
const TrkErrCode&
TrkDifPieceTraj::append(double glen,const TrkDifPieceTraj& other,double& gap)
{
  static TrkErrCode retval;
  retval = TrkErrCode(TrkErrCode::succeed);
// find POCA between the trajectories.
  double mylen = glen;
  HepPoint myend = position(mylen);
  double otherlen = other.lowRange();
  HepPoint otherstart = other.position(otherlen);
  gap = otherstart.distanceTo(myend);
  if(gap > _TOL){ // Don't invoke POCA if the point is too close
    TrkPoca endpoca(other,other.lowRange(),myend,_TOL);
    if(endpoca.status().failure()){
      retval = TrkErrCode(TrkErrCode::fail,7,"poca failure");
      return retval;
    }
    gap = endpoca.doca();
    otherlen = endpoca.flt1();
  }
// get the local trajectory piece at this fltlen
  double loclen(0);
  double piecegap;
// loop over segments in the other trajectory
  for(int itraj = other.trajIndex(otherlen,loclen);
      itraj < other._localtraj.size();itraj++){
// append the piece from the other trajectory
    retval = append(mylen,*(other._localtraj[itraj]),piecegap);
    if(retval.failure())
      break;
// move forward
    mylen = hiRange();
  }
  return retval;
}
// prepend an entire TrkDifTraj
const TrkErrCode&
TrkDifPieceTraj::prepend(double glen,const TrkDifPieceTraj& other,double& gap)
{
  static TrkErrCode retval;
  retval = TrkErrCode(TrkErrCode::succeed);
// find POCA between the trajectories.
  double mylen = glen;
  HepPoint mystart = position(mylen);
  double otherlen = other.hiRange();
  HepPoint otherend = other.position(otherlen);
  gap = otherend.distanceTo(mystart);
  if(gap > _TOL){ // Don't invoke POCA if the point is too close
    TrkPoca endpoca(other,other.hiRange(),mystart,_TOL);
    if(endpoca.status().failure()){
      retval = TrkErrCode(TrkErrCode::fail,7,"poca failure");
      return retval;
    }
    gap = endpoca.doca();
    otherlen = endpoca.flt1();
  }
// get the local trajectory piece at this fltlen
  double loclen(0);
  double piecegap;
// loop over segments in the other trajectory
  for(int itraj = other.trajIndex(otherlen,loclen);
      itraj >= 0;itraj--) {
// prepend the piece from the other trajectory
    retval = prepend(mylen,*(other._localtraj[itraj]),piecegap);
    if(retval.failure())
      break;
// move backwards
    mylen = lowRange();
  }
  return retval;
}
// in-place flight-range inverter.  A traj with range A -> B will end up with
// range -B -> -A, with all the pieces correctly re-ranged as well.
TrkDifPieceTraj&
TrkDifPieceTraj::invert()
{
// make local copies & clear data members
  std::deque<TrkSimpTraj*> trajcopy; trajcopy.swap(_localtraj);
  std::deque<double> rangecopy;      rangecopy.swap(_globalrange);
  assert(_localtraj.size()==0);
  assert(_globalrange.size()==0);
// global range starts at the end of the old range
  _globalrange.push_back(-rangecopy.back());
// loop over the trajectory pieces in backwards order
  for(int itraj=trajcopy.size()-1;itraj>=0;itraj--){
// invert the simptraj and re-append it
    _localtraj.push_back(&(trajcopy[itraj]->invert()));
// set the global range
    _globalrange.push_back(-rangecopy[itraj]);
  }
// reset this traj's range
  double range[2];
  range[0] = -hiRange();
  range[1] = -lowRange();
  Trajectory::setFlightRange(range);
  return *this;
}

// trajectory index from flight range
int
TrkDifPieceTraj::trajIndex(const double& flightdist,double& localflight) const
{
  int index(0); // true when there's only 1 entry
  if (_localtraj.size() > 1) {
    if(validFlightDistance(flightdist)){
//   explicitly check cached value from last call
      if ( _lastIndex >= 0 && _lastIndex <  _localtraj.size()-1 &&
           flightdist > _globalrange[_lastIndex] && flightdist <
           _globalrange[_lastIndex+1] ) {
        index = _lastIndex;
      } else {
//
//  simple binary search algorithm
//
	index = std::upper_bound(_globalrange.begin(), _globalrange.end(), flightdist) - _globalrange.begin() -1;
	if(index>=_localtraj.size()) index=_localtraj.size()-1;
      }
    } else {
//
//  Return the appropriate end if the requested flightdistance is outside the
//  global range
//
      index = (flightdist < lowRange()?0:(_localtraj.size()-1));
    }
  }
  localflight =  localDist(index,flightdist);
// cache value for next time
  _lastIndex = index;
  return index;
}

HepPoint
TrkDifPieceTraj::position(double flightdist) const
{
//
//  First, find the right trajectory piece, then give the local position
//
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(flightdist,localflight);
  return loctraj->position(localflight);
}

Hep3Vector
TrkDifPieceTraj::direction(double flightdist) const
{
//
//  First, find the right trajectory piece, then give the local direction
//
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(flightdist,localflight);
  return loctraj->direction(localflight);
}

double
TrkDifPieceTraj::curvature(double flightdist) const
{
//
//  First, find the right trajectory piece, then give the local curvature
//
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(flightdist,localflight);
  return loctraj->curvature(localflight);
}

Hep3Vector
TrkDifPieceTraj::delDirect(double flightdist) const
{
//
//  First, find the right trajectory piece, then give the local value
//
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(flightdist,localflight);
  return loctraj->delDirect(localflight);
}

void
TrkDifPieceTraj::getInfo(double flightdist,
                         HepPoint& point,
                         Hep3Vector& dir) const
{
//
//  First, find the right trajectory piece, then call the local function
//
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(flightdist,localflight);
  loctraj->getInfo(localflight,point,dir);
}

void
TrkDifPieceTraj::getInfo(double flightdist,
                         HepPoint& point,
                         Hep3Vector& dir,
                         Hep3Vector& deldirect) const
{
//
//  First, find the right trajectory piece, then call the local function
//
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(flightdist,localflight);
  loctraj->getInfo(localflight,point,dir,deldirect);
}

double
TrkDifPieceTraj::distTo1stError(double flightdist,double tol,int dir) const
{
//
//  First, find the right trajectory piece
//
  double localflight(0.0);
  int index = trajIndex(flightdist,localflight);
  const TrkSimpTraj* loctraj = _localtraj[index];
//
//  Ask the local piece for it's dist, and take the minimum of this or the
//  distance to the next trajectory piece
//
  double localdist = loctraj->distTo1stError(localflight,tol,dir);
//
//  Take the minimum of this distance and the distance to the next trajectory
  double dist = localdist;
  if (dir > 0){
    if(index < _localtraj.size()-1)
      dist = std::min(localdist,_globalrange[index+1]-flightdist) + _STEPEPSILON;
  } else {
    if(index > 0)
      dist = std::min(localdist,flightdist - _globalrange[index]) + _STEPEPSILON;
  }
  return dist;
}

double
TrkDifPieceTraj::distTo2ndError(double flightdist,double tol,int dir) const
{
//
//  First, find the right trajectory piece
//
  double localflight(0.0);
  int index = trajIndex(flightdist,localflight);
  const TrkSimpTraj* loctraj = _localtraj[index];
//
//  Ask the local piece for it's dist, and take the minimum of this or the
//  distance to the next trajectory piece
//
  double localdist = loctraj->distTo2ndError(localflight,tol,dir);
//
//  Take the minimum of this distance and the distance to the next trajectory
  double dist = localdist;
  if (dir > 0){
    if(index < _localtraj.size()-1)
      dist = std::min(localdist,_globalrange[index+1]-flightdist) + _STEPEPSILON;
  } else {
    if(index > 0)
      dist = std::min(localdist,flightdist - _globalrange[index]) + _STEPEPSILON;
  }
  return dist;
}


const TrkSimpTraj*
TrkDifPieceTraj::localTrajectory(double flightdist,double& localflight) const
{
//
//  Find and return the right trajectory piece
//
  int index = trajIndex(flightdist,localflight);
  return _localtraj[index];
}
//
//  Re-size the trajectory in the given direction
//
int
TrkDifPieceTraj::resize(double len,trkDirection tdir)
{
  int nremoved(0);
  TrkSimpTraj* trajpiece;
  double oldend,locdist;
  double lrange[2];
  int ipiece = trajIndex(len,locdist);
  switch(tdir) {
  case trkIn:
//
//  Reset the range.  Delete trajectory pieces till we're in range.
//  This will leave the trajectory with a starting range different from
//  0, but it'll still be self-consistent
//
    while(ipiece > 0){
      nremoved++;
      trajpiece = _localtraj.front(); _localtraj.pop_front();
      delete trajpiece;
      oldend = _globalrange.front(); _globalrange.pop_front();
      ipiece = trajIndex(len,locdist);
    }
//  reset the global range
    flightrange[0] = len;
    _globalrange[ipiece] = len;
//  reset the local traj piece range
    lrange[0] = locdist;
    lrange[1] = _localtraj[ipiece]->hiRange();
    _localtraj[ipiece]->setFlightRange(lrange);
    break;
  case trkOut:
    while(ipiece < _localtraj.size() - 1){
      nremoved++;
      trajpiece = _localtraj.back(); _localtraj.pop_back();
      delete trajpiece;
      oldend = _globalrange.back(); _globalrange.pop_back();
      ipiece = trajIndex(len,locdist);
    }
    flightrange[1] = len;
    _globalrange[ipiece+1] = len;
    lrange[0] = _localtraj[ipiece]->lowRange();
    lrange[1] = locdist;
    _localtraj[ipiece]->setFlightRange(lrange);
    break;
  }
  return nremoved;
}
//
//  Print functions
//
void
TrkDifPieceTraj::print(ostream& os) const
{
  os << "TrkDifPieceTraj has " << _localtraj.size() << " pieces "
    << ", total flight range of " << hiRange();
}

void
TrkDifPieceTraj::printAll(ostream& os) const
{
  os << "TrkDifPieceTraj has " << _localtraj.size() << " pieces "
    << ", total flight range from "
     << lowRange() << " to " << hiRange() << endl;
  for(int ipiece=0;ipiece<_localtraj.size();ipiece++){
    TrkSimpTraj* tpiece = _localtraj[ipiece];
    double plow = tpiece->lowRange();
    double phi = tpiece->hiRange();
    os << "Piece " << ipiece << " has global range from "
       << _globalrange[ipiece] << " to " << _globalrange[ipiece+1]
       << " and local range from " <<  plow << " to " << phi << endl;
    HepPoint start = tpiece->position(plow);
    HepPoint end = tpiece->position(phi);
    const HepPoint& refp = tpiece->referencePoint();
    os << "Piece " << ipiece << " starts at point "
       << start.x() <<","<< start.y()<<","<< start.z()
       << " and ends at point "
       << end.x()<<"," << end.y()<<"," << end.z()
       << " and references point "
       << refp.x()<<"," << refp.y()<<"," << refp.z()
       << endl;
  }
}
//
//  Range functions
//
void
TrkDifPieceTraj::setFlightRange(double newrange[2])
{
  double oldend;
  double lrange[2];
  TrkSimpTraj* trajpiece;
//
//  Check for pathological cases
//
  if( newrange[1] > newrange[0] &&
      newrange[0] < flightrange[1] &&
      newrange[1] > flightrange[0] ) {
    resize(newrange[0],trkIn);
    resize(newrange[1],trkOut);
  } else if(newrange[1] < newrange[0] )
    ErrMsg(error) << "cannot resize -- invalid range" << endmsg;
  else {
//
//  Here the new range is completely discontinuous from the
//  old.  We'll just take the first (or last) traj piece,
//  clear out everything else, and reset the ranges.
//
    if( newrange[0] > flightrange[1]){
      while(_localtraj.size()>1){
        trajpiece = _localtraj.front(); _localtraj.pop_front();
        delete trajpiece;
        oldend = _globalrange.front(); _globalrange.pop_front();
      }
    } else {
      while(_localtraj.size()>1){
        trajpiece = _localtraj.back(); _localtraj.pop_back();
        delete trajpiece;
        oldend = _globalrange.back(); _globalrange.pop_back();
      }
    }
    lrange[0] = localDist(0,newrange[0]);
    lrange[1] = localDist(0,newrange[1]);
    flightrange[0] = newrange[0];
    flightrange[1] = newrange[1];
    _globalrange[0] = newrange[0];
    _globalrange[1] = newrange[1];
    _localtraj[0]->setFlightRange(lrange);
  }
}
//
//  KalDeriv functions
//
HepMatrix
TrkDifPieceTraj::derivDeflect(double flightdist,deflectDirection idir) const
{
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(flightdist,localflight);
  return loctraj->derivDeflect(localflight,idir);
}

HepMatrix
TrkDifPieceTraj::derivDisplace(double flightdist,deflectDirection idir) const
{
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(flightdist,localflight);
  return loctraj->derivDisplace(localflight,idir);
}

HepMatrix
TrkDifPieceTraj::derivPFract(double flightdist) const
{
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(flightdist,localflight);
  return loctraj->derivPFract(localflight);
}

const TrkErrCode&
TrkDifPieceTraj::append(TrkSimpTraj* newtraj,double& gap)
{
  static TrkErrCode retval;
  retval = TrkErrCode(TrkErrCode::succeed);
// find out where this piece is relative to the existing pieces.
  HepPoint end = newtraj->position(newtraj->lowRange());
  TrkPoca endpoca(*this,newtraj->lowRange(),end,_TOL);
  if(endpoca.status().success()){
// find the local trajectory pieces at THE end
    double local;
    int index = trajIndex(endpoca.flt1(),local);
// Make sure we're at or beyond the end of the existing trajectory
    if( index == _localtraj.size()-1){
      if(endpoca.flt1() - hiRange() > _SMALLGAP ){
// Fill the gap between the traj pieces with a 0-length copy of the new traj; this
// keeps the bookkeeping straight, provides 'correct' parameters everywhere, 
//  but defines the gap as an invalid range
        double midflt = 0.5*(endpoca.flt1()+lowRange());
        HepPoint mid = position(midflt);
        TrkPoca midpoca(*newtraj,newtraj->lowRange(),mid,_TOL);
        if(midpoca.status().success() && midpoca.flt1() < newtraj->lowRange()){
// create a 0-length copy of the new traj to associate with the midpoint.
// This insures correct coverage of the global flight range
          TrkSimpTraj* gaptraj = (TrkSimpTraj*)(newtraj->clone());
          assert(gaptraj != 0);
          double range[2];
          range[0] = midpoca.flt1();
          range[1] = range[0];
          gaptraj->setFlightRange(range);
// global range of gap traj; from midpoint poca to the newtraj start
          double midrange = newtraj->lowRange() - midpoca.flt1();
// insert the trajectories
          _localtraj.push_back(gaptraj);
          _localtraj.push_back(newtraj);
// append the new global range
          _globalrange.back() = midflt; // override old global fltlen, now goes to midpoint
          _globalrange.push_back(_globalrange.back() + midrange); // gaptraj global range, from midpoint poca to lowRange() of new traj
          _globalrange.push_back(_globalrange.back() + newtraj->range()); // add range of the new traj
// extend the piece-trajs own range
          flightrange[1] = _globalrange.back();
// measure the gap at the midpoint
          gap = midpoca.doca();
        } else
          retval = TrkErrCode(TrkErrCode::fail,17,"poca failure");
      } else {
// just trim the traj to fit
        TrkSimpTraj* endtraj = _localtraj[index];
        static double range[2];
        range[0] = endtraj->lowRange();
        range[1] = local;
        _localtraj.back()->setFlightRange(range); // override local traj range to end at POCA of new piece
        _localtraj.push_back(newtraj);
        _globalrange.back() = endpoca.flt1(); // override old global fltlen, now goes to start of new piece
        _globalrange.push_back(_globalrange.back()+newtraj->range());
        flightrange[1] = _globalrange.back();
        gap = endpoca.doca();
      }
    } else
      retval = TrkErrCode(TrkErrCode::fail,18,"invalid range");
  } else
    retval = TrkErrCode(TrkErrCode::fail,17,"poca failure");
  return retval;
}


const TrkErrCode&
TrkDifPieceTraj::prepend(TrkSimpTraj* newtraj,double& gap)
{
  static TrkErrCode retval;
  retval = TrkErrCode(TrkErrCode::succeed);
// find out where this piece is relative to the existing pieces.
  HepPoint end = newtraj->position(newtraj->hiRange());
  TrkPoca poca(*this,newtraj->hiRange(),end,_TOL);
  if(poca.status().success()){
// find the local trajectory pieces at either end
    double local;
    int index = trajIndex(poca.flt1(),local);
    TrkSimpTraj* oldtraj = _localtraj[index];
// Make sure we're beyond the end of the existing trajectory
    if( index == 0 && poca.flt1() < lowRange()){
// we want to split the gap between the traj pieces evenly,
      double gapstart = poca.flt1();
      double gapend = globalDist(index,oldtraj->hiRange());
      double gapmid = (gapend+gapstart)/2.0;
      HepPoint mid = position(gapmid);
// approximate initial (local) fltlen for poca
      double locmid = newtraj->hiRange() - gapstart + gapmid;
      TrkPoca midpoca(*newtraj,locmid,mid,_TOL);
      if(midpoca.status().success()){
// create a 0-length copy of the new traj to associate with the midpoint.
// This insures correct coverage of the global flight range
        TrkSimpTraj* gaptraj = (TrkSimpTraj*)(newtraj->clone());
        assert(gaptraj != 0);
        double range[2];
        range[0] = midpoca.flt1();
        range[1] = range[0];
        gaptraj->setFlightRange(range);
// insert the trajectories
        _localtraj.push_front(gaptraj);
        _localtraj.push_front(newtraj);
// append the new global range
        _globalrange.push_front(gapmid);
        _globalrange.push_front(gapstart-newtraj->range());
// extend the piece-trajs own range
        flightrange[0] = _globalrange.front();
// measure the gap at the midpoint
        gap = mid.distanceTo(newtraj->position(midpoca.flt1()));
      } else
        retval = TrkErrCode(TrkErrCode::fail,17,"poca failure");
    } else
      retval = TrkErrCode(TrkErrCode::fail,18,"invalid range");
  } else
    retval = TrkErrCode(TrkErrCode::fail,17,"poca failure");
  return retval;
}

bool
TrkDifPieceTraj::locallyValid(double glen,double tol) const
{
  double localflight(0.0);
  const TrkSimpTraj* loctraj = localTrajectory(glen,localflight);
  return loctraj != 0 &&
    localflight-tol <= loctraj->hiRange()  &
    localflight+tol >= loctraj->lowRange();
}

bool
TrkDifPieceTraj::operator == (const TrkDifPieceTraj& other) const {
  bool retval(true);
  if(retval) retval = _globalrange == other._globalrange;
  if(retval) retval = _localtraj.size() == other._localtraj.size();
// loop over component trajs
  std::deque<TrkSimpTraj*>::const_iterator miter = _localtraj.begin();
  std::deque<TrkSimpTraj*>::const_iterator oiter = other._localtraj.begin();
  while(retval && miter != _localtraj.end() && oiter != other._localtraj.end()){
    retval = (*oiter)->parameters()->parameter() == (*miter)->parameters()->parameter() &&
      (*oiter)->parameters()->covariance() == (*miter)->parameters()->covariance();
    miter++;oiter++;
  }
  return retval;
}
