// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalPairSite.cc,v 1.19 2004/08/06 06:12:48 bartoldu Exp $
//
//    class KalPairSite.  
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/21/98
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include "KalmanTrack/KalPairSite.hh"
#include "KalmanTrack/KalPairVisitor.hh"
#include "FrameLogger/AppErrLog.hh"
#include "BbrGeom/BbrVectorErr.hh"
#include "BbrGeom/BbrDoubleErr.hh"
#include "ErrLogger/ErrLog.hh"
using std::ostream;

KalPairSite::KalPairSite(const TrkDifPieceTraj* ptraj,
			 const TrkParams& constraint, bool* constrainParams,
			 double fltlen, KalPairConduit* conduit) : 
  KalConstraint(ptraj, constraint, constrainParams, fltlen, KalSite::pairSite),
  _conduit(conduit)
{
}
KalPairSite::~KalPairSite(){;}

KalPairSite*
KalPairSite::clone(const KalRep* krep) const {
  // Don't want to be clonable.  The conduit interface would be messed up
  return 0;
}

bool
KalPairSite::update(const TrkDifPieceTraj* newtraj, double newmom) {

  // See if the conduit has some new information for this site
  assert(0 != _conduit);
  if (_conduit->checkLAM(this)) {
    double fltlen = _conduit->getFltLen(this);
    changeConstraint(_conduit->getConstraint(this));
    //    _fltchi2 = _conduit->getFltChi2(this);
    setTraj(newtraj, fltlen);
    reset();
  }
  // Call down to KalConstraint::update
  return KalConstraint::update(newtraj, newmom);
}

bool
KalPairSite::process(const KalSite* prevsite, trkDirection idir) {

  assert (0 != _conduit);
  // Send the inward fit results to the conduit
  if (idir == trkIn) {
    _conduit->uploadParams(prevsite->filterParameters(idir), this);
  }
  // If the constraint status is non-zero, this hasn't been updated
  // yet (happens on first call).  So, just return without actually
  // doing anything.
  return KalConstraint::process(prevsite, idir);
}
void
KalPairSite::printAll(ostream& os) const {
   os << " KalPairSite ";
   KalConstraint::printAll(os);
}

void
KalPairSite::changeConstraint(const KalParams& params) {

  _cparams = params;
  HepSymMatrix cmat = params.covarianceMatrix();
  int ierr(0);
  cmat.invert(ierr);
  if(ierr == 0) {
// mask off the information for un-constrained parameters
    unsigned npar = _cparams.nPar();
    for(unsigned ipar=0;ipar<npar;ipar++){
      for(unsigned jpar=0;jpar<=ipar;jpar++){
	if((!isConstrained(ipar)) || (!isConstrained(jpar)) )
	  cmat.fast(ipar+1,jpar+1) = 0.0;
      }
    }
    HepVector cvec(cmat*params.parameterVector());
    _cweight = KalWeight(cvec,cmat);
  }

  return;
}
