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
//	20020522  M. Kelsey -- Remove assert() from resid(HOT*...).  Replace
//		  with sanity checks on HOT/Rep association and whether HOT
//		  has already-computed residual.  Return value used to
//		  flag sanity checks and "trustworthiness" of results.
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "BaBar/Constants.hh"
#include "BField/BFieldFixed.hh"
#include <assert.h>
#include <algorithm>
#include <iostream>
#include "TrkBase/TrkRep.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkHotListFull.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkFunctors.hh"
#include "TrkBase/TrkErrCode.hh"
#include "difAlgebra/DifPoint.hh"
#include "difAlgebra/DifVector.hh"
#include "difAlgebra/DifIndepPar.hh"
#include "ProbTools/ChisqConsistency.hh"
#include "ErrLogger/ErrLog.hh"
#include "TrkBase/TrkExchangePar.hh"
using std::cout;
using std::endl;

TrkRep::TrkRep(TrkParticle const& hypo,bool createHotList)
  : _tpart(hypo),
    _betainv(-999999.),
    _hotList( createHotList?new TrkHotListFull:0 )
{
}

TrkRep::TrkRep(const TrkHotList& hotlist, 
               TrkParticle const& hypo)
  : _tpart(hypo),
    _betainv(-999999.),
    _hotList( hotlist.clone(TrkBase::Functors::cloneHot(this)) )
{
}

TrkRep::TrkRep(TrkHotList& hotlist, 
               TrkParticle const& hypo, bool stealHots)
  : _tpart(hypo),
    _betainv(-999999.),
    _hotList( stealHots? new TrkHotListFull(hotlist,setParent(this))
                       : hotlist.clone(TrkBase::Functors::cloneHot(this)) )
{
}

TrkRep::TrkRep(const TrkHotList* hotlist, 
               TrkParticle const& hypo)
  : _tpart(hypo),
    _betainv(-999999.),
    _hotList( hotlist!=0?
                  hotlist->clone(TrkBase::Functors::cloneHot(this)):
                  new TrkHotListFull )
{
}

TrkRep::TrkRep(TrkHotList* hotlist, 
               TrkParticle const& hypo,bool takeownership)
  : _tpart(hypo),
    _betainv(-999999.)
{
  if (!takeownership) {
    _hotList.reset( hotlist!=0?
                    hotlist->clone(TrkBase::Functors::cloneHot(this)):
                    new TrkHotListFull );
  } else {
    assert(hotlist!=0);
    _hotList.reset( hotlist->resetParent(setParent(this)) );
  }
}

// copy ctor
TrkRep::TrkRep(const TrkRep& oldRep,  TrkParticle const& hypo) :
  TrkFitStatus(oldRep),
    _tpart(hypo),
    _betainv(-999999.)
{
  // Hots and hotlist have to be cloned in the derived classes
}

TrkRep&
TrkRep::operator= (const TrkRep& right)
{
  if(&right != this){
    _tpart=right._tpart;
    _betainv=right._betainv;
    _hotList.reset( right._hotList->clone(this) );
    TrkFitStatus::operator=(right);
  }
  return *this;
}

TrkRep::~TrkRep()
{
}

bool
TrkRep::operator== (const TrkRep& rhs)
{
  return (&rhs == this);
}

void
TrkRep::addHot(TrkHitOnTrk *newHot)
{
  newHot->setParent(this);
  if (newHot->isActive()) setCurrent(false);
  hotList()->append(newHot);
}

void
TrkRep::removeHot(TrkHitOnTrk *theHot)
{
  if(theHot->isActive()) setCurrent(false);     // fit no longer current
  hotList()->remove(theHot);
}

void
TrkRep::activateHot(TrkHitOnTrk* hot)
{
  if(!hot->isActive()){
// make sure this is my hot we're talking about
    if(this == hot->getParentRep()){
      setCurrent(false);
// actually activate the hot; this is now the rep's job
      hot->setActive(true);
    }
  }
}

void
TrkRep::deactivateHot(TrkHitOnTrk* hot)
{
  if(hot->isActive()){
// make sure this is my hot we're talking about
    if(this == hot->getParentRep()){
      setCurrent(false);
// actually deactivate the hot; this is now the rep's job
      hot->setActive(false);
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
  static double cinv = 1./Constants::c;
  // Initialize cache
  if (_betainv < 0.0) {
    double mass2 = particleType().mass();
    mass2 = mass2 * mass2;
    double ptot2 = momentum(0.).mag2();
    assert(ptot2 != 0.0);
    _betainv = sqrt( (ptot2 +  mass2)/ ptot2);
  }
  double tof = fltL * _betainv * cinv;
  return _trkt0.t0() + tof;
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
  return hotList()->startFoundRange();
}

double
TrkRep::endFoundRange() const
{
  return hotList()->endFoundRange();
}

void
TrkRep::updateHots()
{
  setCurrent(false);
  hotList()->updateHots();
}

int
TrkRep::nActive() const
{
  return hotList()->nActive();
}

bool
TrkRep::resid(const TrkHitOnTrk *h,
              double& residual, double& residErr,
              bool exclude) const
{
  assert (h != 0);
  if (h->parentRep() != this) return false;	// HOT must belong to Rep
  if (!h->hasResidual()) return false;		// Residual must be available
  if (exclude) return false;  			// FIXME: Can't do unbiased residuals (yet!)

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

BField const&
TrkRep::bField() {
// this should be built from the Mu2e field: FIXME!!!
  static BField* bfield = new BFieldFixed(0.0,0.0,1.0);
  return *bfield;
}
