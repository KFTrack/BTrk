//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPocaBase.cc,v 1.30 2006/03/25 15:15:55 brownd Exp $
//
// Description:
//     
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
// Modifications:
// - Jan 2003 (WDH) 
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkPocaBase.hh"
#include "TrkBase/TrkPoca.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/HepPoint.h"
#include "BbrGeom/Trajectory.hh"
#include <math.h>
#include <assert.h>
#include "ErrLogger/ErrLog.hh"

#include <iomanip>
#include <algorithm>

TrkPocaBase::TrkPocaBase(double f1, double f2, double prec)
        :  _precision(prec), _flt1(f1), _flt2(f2),_status(TrkErrCode::fail)
{
  assert(prec > 0.);
}

TrkPocaBase::TrkPocaBase()
        :  _precision(1.0e-5), _flt1(0.0), _flt2(0.0),_status(TrkErrCode::fail)
{
}

TrkPocaBase::TrkPocaBase(const TrkPocaBase& other)
  :  _precision(other._precision), _flt1(other._flt1), _flt2(other._flt2),
     _status(other._status)
{
}

TrkPocaBase&
TrkPocaBase::operator = (const TrkPocaBase& other) {
  if(this != &other){
    _precision = other._precision;
    _flt1 = other._flt1;
    _flt2 = other._flt2;
    _status = other._status;
  }
  return *this;
}

void 
TrkPocaBase::minimize(TrkPocaTraj& ptraj1,
                      TrkPocaTraj& ptraj2)
{
  // Last revision: Jan 2003, WDH
  const int maxnOscillStep = 5 ;
  const int maxnDivergingStep = 5 ;
  const int maxnStuck = 3 ;

  // initialize
  _status = TrkErrCode::succeed;
  
  static HepPoint newPos1, newPos2 ;
  double delta(0), prevdelta(0) ;
  int nOscillStep(0) ;
  int nDivergingStep(0) ;
  int nStuck(0);
  bool finished = false ;
  int istep(0) ;

  for (istep=0; istep < _maxTry && !finished; ++istep) {
    double prevflt1  = ptraj1._flt;
    double prevflt2  = ptraj2._flt ;
    double prevprevdelta = prevdelta ;
    prevdelta = delta;
    
    stepTowardPoca(ptraj1, ptraj2);
    if( status().failure() ) {
      // failure in stepTowardPoca
      finished=true ;
    } else {
      newPos1 = ptraj1._traj.position(ptraj1._flt);
      newPos2 = ptraj2._traj.position(ptraj2._flt);
      delta = (newPos1 - newPos2).mag();
      double step1 = ptraj1._flt - prevflt1;
      double step2 = ptraj2._flt - prevflt2;
      int pathDir1 = (step1 > 0.) ? 1 : -1;
      int pathDir2 = (step2 > 0.) ? 1 : -1;
      
      // Can we stop stepping?
      double distToErr1 = ptraj1._traj.distTo1stError(prevflt1, precision(), pathDir1);
      double distToErr2 = ptraj2._traj.distTo1stError(prevflt2, precision(), pathDir2);
      
      // converged if very small steps, or if parallel
      finished = 
	(fabs(step1) < distToErr1 && fabs(step2) < distToErr2 ) ||
	(status().success() == 3) ;
 
      // we have to catch some problematic cases
      if( !finished && istep>2 && delta > prevdelta) {
      // we can get stuck if a flt range is restricted
	if(ptraj1._rflt && step1==0.0 ||
	   ptraj2._rflt && step2==0.0) {
	  if(++nStuck > maxnStuck){
	  // downgrade to a point poca
	    if(ptraj2._rflt)
	      minimize(ptraj1,newPos2);
	    else
	      minimize(ptraj2,newPos1);
	    _status.setSuccess(22,"Stuck poca.");
	    finished = true;
	  }
	} else if( prevdelta > prevprevdelta) {
	  // diverging
	  if(++nDivergingStep>maxnDivergingStep) { 
	    _status.setFailure(2) ; // code for `Failed to converge'
	    finished = true ;
	  }
	} else {
 	  nDivergingStep=0;
	  // oscillating
	  if(++nOscillStep>maxnOscillStep) {
	    // bail out of oscillation. since the previous step was
	    // better, use that one.
	    ptraj1._flt = prevflt1 ;
	    ptraj2._flt = prevflt2 ;
	    _status.setSuccess(21, "Oscillating poca.") ;
	    finished = true ;
	  } else {
	    // we might be oscillating, but we could also just have
	    // stepped over the minimum. choose a solution `in
	    // between'.
	    setFlt(prevflt1 + 0.5*step1,ptraj1);
	    setFlt(prevflt2 + 0.5*step2,ptraj2);
	    newPos1 = ptraj1._traj.position(ptraj1._flt) ;
	    newPos2 = ptraj2._traj.position(ptraj2._flt) ;
	    delta = (newPos1 - newPos2).mag() ;
	  }
	}
      } 
    }
  }
  if(!finished) _status.setSuccess(2) ; // code for 'not converged' (yet)
}

TrkPocaBase::TrkPocaBase(double f1, double prec) 
  : _precision(prec), _flt1(f1), _flt2(0), _status(TrkErrCode::fail)
{
}

void
TrkPocaBase::minimize(TrkPocaTraj& ptraj,const HepPoint& pt )
{
  _status=TrkErrCode::succeed;
  int pathDir = 1;  // which way are we stepping (+/- 1)

  int nTinyStep = 0;  // number of consecutive tiny steps -- oscillation test
  int nOscills = 0;
  double fltLast = 0., fltBeforeLast = 0.; // another check for oscillation
  for (int i = 0; i < _maxTry; i++) {
    fltLast = ptraj._flt;
    stepToPointPoca(ptraj, pt);
    if (status().failure()) return;
    double step = ptraj._flt - fltLast;
    pathDir = (step > 0.) ? 1 : -1;
    // Can we stop stepping?
    double distToErr = ptraj._traj.distTo1stError(fltLast, precision(), pathDir);
    bool mustStep = (fabs(step) > distToErr && step != 0.);
    // Crude test for oscillation (around cusp point of piecewise traj, I hope)
    if (fabs(step) < 0.5 * precision()) {
      nTinyStep++;
    } else {
      nTinyStep = 0;
      if (i > 1) {
        if (fabs(step) >= fabs(fltBeforeLast-fltLast) &&
            fabs(fltBeforeLast-ptraj._flt) <= fabs(step)) {
          nOscills++;
          double halfway = (fltBeforeLast + fltLast) / 2.;
          if ((ptraj._traj.position(ptraj._flt) - pt).mag() >
              (ptraj._traj.position(halfway) - pt).mag()) ptraj._flt = halfway;
        }
      }
    }
    if (nTinyStep > 3) mustStep = false;
    if (nOscills > 2) {
      mustStep = false;
      ErrMsg(warning) << "Alleged oscillation detected. "
                      << step << "  " << fltLast-fltBeforeLast
                      << "  " << i << "  " << endmsg;
    }
    if (!mustStep) return;
    fltBeforeLast = fltLast;
  }
  // Ran off the end of the loop
  _status.setFailure(2);
}


TrkPocaBase::~TrkPocaBase() 
{}


void 
TrkPocaBase::stepTowardPoca(TrkPocaTraj& ptraj1,
			    TrkPocaTraj& ptraj2)
{ 
  // Last revision: Jan 2003, WDH
  
  // A bunch of unsightly uninitialized variables:
  static Hep3Vector dir1, dir2;
  static Hep3Vector delDir1, delDir2;
  static HepPoint pos1, pos2;

  ptraj1._traj.getInfo(ptraj1._flt, pos1, dir1, delDir1);
  ptraj2._traj.getInfo(ptraj2._flt, pos2, dir2, delDir2);
  Hep3Vector delta = pos1 - pos2;
  double ua = -delta.dot(dir1);
  double ub =  delta.dot(dir2);
  double caa = dir1.mag2() + delta.dot(delDir1);
  double cbb = dir2.mag2() - delta.dot(delDir2);
  double cab = -dir1.dot(dir2);
  double det = caa * cbb - cab * cab;
  
  if(det<0) {
  // get rid of second order terms
    caa = dir1.mag2() ;
    cbb = dir2.mag2() ;
    det = caa * cbb - cab * cab;
  }
  
  if ( det < 1.e-8) {
    // If they are parallel (in quadratic approximation) give up
    _status.setSuccess(3);
    return;
  }
  
  double df1 = (ua * cbb - ub * cab)/det;
  int pathDir1 = (df1 > 0) ? 1 : -1;
  double df2 = (ub * caa - ua * cab)/det;
  int pathDir2 = (df2 > 0) ? 1 : -1;
  
  // Don't try going farther than worst parabolic approximation will
  // allow: Since ` _extrapToler' is large, this cut effectively only
  // takes care that we don't make large jumps past the kink in a
  // piecewise trajectory.

  double distToErr1 = ptraj1._traj.distTo2ndError(ptraj1._flt, _extrapToler, pathDir1);
  double distToErr2 = ptraj2._traj.distTo2ndError(ptraj2._flt, _extrapToler, pathDir2);
  
  // Factor to push just over border of piecewise traj (essential!)
  const double smudge = 1.01 ; 
  if( fabs(df1) > smudge*distToErr1 ) {
    // choose solution for which df1 steps just over border
    df1 = smudge*distToErr1 * pathDir1 ;
    // now recalculate df2, given df1:
    df2 = (ub - df1*cab)/cbb ;
  }

  if( fabs(df2) > smudge*distToErr2 ) {
    // choose solution for which df2 steps just over border
    df2 = smudge*distToErr2 * pathDir2 ;
    // now recalculate df1, given df2:
    df1 = (ua - df2*cab)/cbb ;
    // if still not okay,
    if( fabs(df1) > smudge*distToErr1 ) {
      df1 = smudge*distToErr1 * pathDir1 ;
    }
  }
  
  setFlt(ptraj1._flt+df1,ptraj1);
  setFlt(ptraj2._flt+df2,ptraj2);

  // another check for parallel trajectories
  if (fabs(ptraj1._flt) > _maxDist && fabs(ptraj2._flt) > _maxDist) 
    _status.setSuccess(3) ;
}

void 
TrkPocaBase::stepToPointPoca(TrkPocaTraj& ptraj, const HepPoint& pt) 
{
// Unsightly uninitialized variables:
  static Hep3Vector dir, delDir;
  static HepPoint trajPos;

  ptraj._traj.getInfo(ptraj._flt, trajPos, dir, delDir);
  Hep3Vector delta = trajPos - pt;
  double denom = 1. + delta.dot(delDir);
  if (fabs(denom)*_maxDist < 1. ) {
    _status.setFailure(11, "TrkPoca::ambiguous tight looper.");
    return;
  }
  double df = -delta.dot(dir) / fabs(denom);
  int pathDir = (df > 0.) ? 1 : -1;

  // Don't try going farther than worst parabolic approximation will allow:
  double distToErr = ptraj._traj.distTo2ndError(ptraj._flt, _extrapToler, pathDir);
  if (fabs(df)>distToErr) df = (df>0?distToErr:-distToErr);
  // Make the step slightly longer -- prevents quitting at kinks
  df += 0.001 * pathDir * precision();
  setFlt(ptraj._flt+df,ptraj);
}

double TrkPocaBase::_maxDist = 1.e7;

int TrkPocaBase::_maxTry = 500;

double TrkPocaBase::_extrapToler = 2.;

void
TrkPocaBase::setFlt(double flt,TrkPocaTraj& ptraj) {
  if(!ptraj._rflt)
    ptraj._flt = flt;
  else
    ptraj._flt = std::max(std::min(ptraj._traj.hiRange(),flt),ptraj._traj.lowRange());
}
