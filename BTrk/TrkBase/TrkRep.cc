//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkRep.cc,v 1.73 2004/11/04 22:21:51 raven Exp $
//
// Description:
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//
// Revision History (started 2002/05/22)
//	20020522  M. Kelsey -- Remove assert() from resid(TrkHit*...).  Replace
//		  with sanity checks on TrkHit/Rep association and whether TrkHit
//		  has already-computed residual.  Return value used to
//		  flag sanity checks and "trustworthiness" of results.
//------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/Constants.hh"
#include "BTrk/BField/BFieldFixed.hh"
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/TrkBase/TrkDifTraj.hh"
#include "BTrk/TrkBase/TrkHit.hh"
#include "BTrk/TrkBase/TrkFunctors.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/difAlgebra/DifPoint.hh"
#include "BTrk/difAlgebra/DifVector.hh"
#include "BTrk/difAlgebra/DifIndepPar.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BaBar/ErrLog.hh"
#include "BTrk/BbrGeom/BbrPointErr.hh"
#include "BTrk/TrkBase/HelixParams.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Matrix.h"

#include <assert.h>
#include <algorithm>
#include <iostream>

using std::cout;
using std::endl;
using namespace CLHEP;

TrkRep::TrkRep(const TrkHitVector& hitlist, TrkParticle const& hypo)
  : _tpart(hypo), _hitvec( hitlist)
{
  sortHits();
}

// cleanup TrkHits
TrkRep::~TrkRep()
{
  for(auto ihit=_hitvec.begin();ihit!=_hitvec.end();++ihit)
    delete *ihit;
}

void
TrkRep::sortHits() {
// sort hits by flightlength
 std::sort(_hitvec.begin(),_hitvec.end(),hitsort());
}

// disallow
TrkRep&
TrkRep::operator = (const TrkRep& rhs) {
  assert(false);
  return *this;
}

bool
TrkRep::operator== (const TrkRep& rhs)
{
  return (&rhs == this);
}

void
TrkRep::addHit(TrkHit* newTrkHit)
{
  newTrkHit->setParent(this);
  if (newTrkHit->isActive()) setCurrent(false);
  _hitvec.push_back(newTrkHit);
  sortHits();
}

void
TrkRep::removeHit(TrkHit* theTrkHit)
{
  if(theTrkHit->isActive()) setCurrent(false);     // fit no longer current
  auto ifnd = std::find(_hitvec.begin(),_hitvec.end(),theTrkHit);
  if(ifnd != _hitvec.end())
    _hitvec.erase(ifnd);

}

void
TrkRep::activateHit(TrkHit* thit)
{
  if(!thit->isActive()){
// make sure this is my thit we're talking about
    if(this == thit->getParentRep()){
      setCurrent(false);
// actually activate the thit; this is now the rep's job
      thit->setActive(true);
    }
  }
}

void
TrkRep::deactivateHit(TrkHit* thit)
{
  if(thit->isActive()){
// make sure this is my thit we're talking about
    if(this == thit->getParentRep()){
      setCurrent(false);
// actually deactivate the thit; this is now the rep's job
      thit->setActive(false);
    }
  }
}

HepPoint
TrkRep::position(double fltL) const
{
  return traj().position(fltL);
}

Hep3Vector
TrkRep::direction(double fltL) const
{
  return traj().direction(fltL);
}

double
TrkRep::arrivalTime(double fltL) const
{
// average momentum between flt0 and this flight.  Should really integrate, FIXME!!!
  double mom0 = momentum(_flt0).mag();
  double momend = momentum(fltL).mag();
  double avgmom = 0.5*(mom0+momend);
  double beta = _tpart.beta(avgmom);
  return _trkt0.t0() + (fltL-_flt0)/(beta*Constants::c);
}

BbrPointErr
TrkRep::positionErr(double fltL) const
{
  static DifPoint posD;
  static DifVector dirD;
  traj().getDFInfo2(fltL, posD, dirD);
  HepMatrix err = posD.errorMatrix( posD.x.indepPar()->covariance() );
  HepPoint point(posD.x.number(), posD.y.number(), posD.z.number());
  BbrError symErr(3);
  symErr.assign(err);

  //  if (ErrLogging(debugging)) {
  if (false) {
    ErrMsg(routine) << "Pos " << err.num_row() << " " << err.num_col()
                    << endl
                    << "output:" << endl

      //    << err(1,1) << endl
      //    << err(2,1) << "  " << err(2,2) << endl
      //    << err(3,1) << "  " << err(3,2) << "  " << err(3,3) << endl
                    << "x deriv: " << endl
                    << posD.x.derivatives() << endl
                    << "y deriv: " << endl
                    << posD.y.derivatives() << endl
                    << endmsg;
    //  }

    Hep3Vector pointDir(point.x(), point.y());
    double dirMag = pointDir.mag();
    double dist = 5.e-3;
    double delx = dist * point.x() / dirMag;
    double dely = dist * point.y() / dirMag;
    int ierr = 0;
    HepMatrix weight = err.inverse(ierr);
    double chisq =     weight(1,1) * delx * delx +
      2 * weight(2,1) * delx * dely +
      weight(2,2) * dely * dely;
    cout << point << endl;
    cout << symErr << endl;
    cout << "delta: " << delx << "  " << dely << endl;
    cout << "chisq: " << chisq << endl;
    double phi0 = helix(fltL).phi0();
    delx = dist * cos(phi0);
    dely = dist * sin(phi0);
    chisq =            weight(1,1) * delx * delx +
                   2 * weight(2,1) * delx * dely +
                       weight(2,2) * dely * dely;
    cout << "delta: " << delx << "  " << dely << endl;
    cout << "chisq: " << chisq << endl;

    cout << endl << endl;
  }
  return BbrPointErr(point, symErr);
}

BbrVectorErr
TrkRep::directionErr(double fltL) const
{
  static DifPoint posD;
  static DifVector dirD;
  traj().getDFInfo2(fltL, posD, dirD);
  BbrError symErr(3);
  symErr.assign( dirD.errorMatrix( dirD.x.indepPar()->covariance() ));
  Hep3Vector dir(dirD.x.number(), dirD.y.number(), dirD.z.number());
  return BbrVectorErr(dir, symErr);
}

double
TrkRep::startValidRange() const
{
  return traj().lowRange();
}

double
TrkRep::endValidRange() const
{
  return traj().hiRange();
}

double
TrkRep::startFoundRange() const
{
  return _hitvec.front()->fltLen();
}

double
TrkRep::endFoundRange() const
{
  return _hitvec.back()->fltLen();
}

int
TrkRep::nActive() const
{
  int retval(0);
  for(auto ihit=_hitvec.begin();ihit!=_hitvec.end();++ihit){
    if((*ihit)->isActive())++retval;
  }
  return retval;
}

int
TrkRep::nHits() const {
  return _hitvec.size();
}

bool
TrkRep::resid(const TrkHit* h,
              double& residual, double& residErr,
              bool exclude) const
{
  assert (h != 0);
  if (h->parentRep() != this) return false;	// TrkHit must belong to Rep
  if (!h->hasResidual()) return false;		// Residual must be available
  if (exclude) return false;  			// Can't do unbiased residual in base class

  residual=h->residual();
  residErr=h->hitRms();
  return true;
}

ChisqConsistency
TrkRep::chisqConsistency() const {
  if(fitValid())
    return ChisqConsistency(chisq(),nDof());
  else
    return ChisqConsistency();
}

