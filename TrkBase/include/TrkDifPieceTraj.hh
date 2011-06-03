// ----------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkDifPieceTraj.hh,v 1.35 2004/08/06 06:31:39 bartoldu Exp $
//
//  Description:
//  Differential form of the Piecewise compound trajectory class.  Note that this
//  is now a FULLY IMPLEMENTED CLASS, which completely and generally satisfies
//  the Kalman track interface.  Note also that it's now possible to have a
//  heterogenous collection of TrkSimpTrajs (IE you can append difLines on to
//  difHelixes).
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/17/97
//-----------------------------------------------------------------------------
#ifndef TRKDIFPIECETRAJ_HH
#define TRKDIFPIECETRAJ_HH
#include "TrkBase/TrkSimpTraj.hh"
#include "TrkBase/TrkKalDeriv.hh"
#include "TrkBase/TrkDirection.hh"
#include <vector>
#include <deque>

#include <iosfwd>
class TrkErrCode;

class TrkDifPieceTraj :  public TrkDifTraj , public TrkKalDeriv {
public:

//***
//  Constructors
//***
  TrkDifPieceTraj(const TrkSimpTraj&,const double lowlim ,const double hilim );
// same as above, except it _TAKES OWNERSHIP_ of the traj
  TrkDifPieceTraj(TrkSimpTraj*,const double lowlim ,const double hilim );
  TrkDifPieceTraj(const TrkDifPieceTraj& ); // copy constructor
// construct from a vector of simptrajs (which will be cloned).  The first segment will define
// the global flightlength origin, and subsequent global lengths will be computed using POCA.
  TrkDifPieceTraj(const std::vector<TrkSimpTraj*>& trajs);
// construct from a vector of simptrajs, which are TAKEN FOR OWNERSHIP.  The first segment will define
// the global flightlength origin, and subsequent global lengths will be computed using POCA.
  TrkDifPieceTraj(std::vector<TrkSimpTraj*>& trajs,double gstart=0.0);  
  virtual ~TrkDifPieceTraj();
  TrkDifPieceTraj& operator =(const TrkDifPieceTraj&);
// validation test
  bool operator == (const TrkDifPieceTraj& other) const;
// inversion function, which keeps the same 'trajectory' but reverses the flightlength.
// this is done in-place, and returns itself
  TrkDifPieceTraj& invert();
  TrkDifPieceTraj* clone() const  {return new TrkDifPieceTraj(*this); }

//***
//  Differential functions: all are implemented here as call downs
//***
  void getDFInfo(double fltLen, DifPoint& pos, DifVector& direction,
		 DifVector& delDirect) const;
  void getDFInfo2(double fltlen, DifPoint& pos, DifVector& direction) const;

//***
// Normal trajectory functions: these are all fully implemented as call downs
//***
  HepPoint position( double ) const;
  Hep3Vector direction( double ) const;
  double curvature( double f = 0. ) const;
  Hep3Vector delDirect( double ) const;
  void  getInfo(double fltLen, HepPoint& , Hep3Vector& direction) const;
  void  getInfo(double fltLen, HepPoint& , Hep3Vector& direction,
                Hep3Vector& delDirect) const;
  double distTo1stError(double s, double tol, int pathDir) const;
  double distTo2ndError(double s, double tol, int pathDir) const;

//***
//  DifPieceTraj specific functions
//***
  const TrkSimpTraj* localTrajectory(double,double&) const; // const version
// test validity of the local range given a global flightlength.
  bool locallyValid(double glen,double tol=0.0) const;
//  Overwrite the range functions, to allow local ranges to be
//  updated as well

  void setFlightRange(double newrange[2]);

//***
//  print functions
//***

  void print(std::ostream& os) const;
  void printAll(std::ostream& os) const;

//  append/prepend new trajectory pieces.  These will adjust the local range of the
//  appended trajectory to match as best as possible (smallest gap) with the existing
//  trajectory, at the specified global flight length.  The return value is the magnitude
//  of the gap at the append point.  Note that this function will remove existing trajectories
//  from the piecetraj if those are found to overlap with the specified global range.
//
  const TrkErrCode& append(double gfltlen,const TrkSimpTraj&,double& gap);
  const TrkErrCode& prepend(double gfltlen,const TrkSimpTraj&,double& gap);
// ownership-taking versions of the above
  const TrkErrCode& append(double gfltlen,TrkSimpTraj*,double& gap);
  const TrkErrCode& prepend(double gfltlen,TrkSimpTraj*,double& gap);
// append entire piecetrajs together.  This will adjust the range of the first piece as
// with the single piece append/prepend, but will leave the rest of the piecetraj as-is.
  const TrkErrCode& append(double,const TrkDifPieceTraj&,double& gap);
  const TrkErrCode& prepend(double,const TrkDifPieceTraj&,double& gap);
//  Insert the traj 'as-is', without trying to adjust the
//  local trajectories flight length to maximize continuity.  The traj will be inserted
//  at the global range closest to its match with the nearest trajectory neighbors.
//  The input trajectory should be disjoint from the existing pieces, otherwise existing
//  pieces may get deleted.  The input trajectory must have a valid flight range.  'Empty
//  spaces' in the global flightrange will be filed by clones of existing trajectories.
//  Index and localTraj functions will operate on the resultant piecetraj to find the
//  closest piece to the specified flightlength, although the returned local flight length
//  may be outside the localtrajs range.
//  This function _TAKES OWNERSHIP_ of the input traj.
//
  const TrkErrCode& append(TrkSimpTraj* traj,double& gap);
  const TrkErrCode& prepend(TrkSimpTraj* traj,double& gap);
//***
//  Real versions of the base class derrivative functions
//***
  HepMatrix derivDeflect(double fltlen,deflectDirection) const;
  HepMatrix derivDisplace(double fltlen,deflectDirection idir) const;

  HepMatrix derivPFract(double fltlen) const;
protected:
  int resize(double len,trkDirection); // extend/shorten the trajectory.
// find the trajectory index given the global flight distance:
  int trajIndex(const double& global,double& local) const;
// global to local translation:
  double localDist(int index,double globdist) const {
    return _localtraj[index]->lowRange() + globdist - _globalrange[index];
  }
// local to global translation
  double globalDist(int index,double locdist) const {
    return _globalrange[index] + locdist -  _localtraj[index]->lowRange();
  }
// Data Members
  std::deque<TrkSimpTraj*> _localtraj; // Ordered vector of simptrajs
  std::deque<double> _globalrange; // local/global flight range conversion
  mutable int _lastIndex;  // Cache of last index used to access traj -- speeds up
                   // searching (iff the same piece is queried repeatedly)
};
#endif
