//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalPairConduit.cc,v 1.10 2004/06/21 22:07:53 brownd Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Doug Roberts
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "AbsEnv/AbsEnv.hh"
#include "BbrGeom/BbrDoubleErr.hh"
#include "BField/BField.hh"
#include "KalmanTrack/KalPairConduit.hh"
#include "KalmanTrack/KalPairVisitor.hh"
#include "TrkBase/TrkDifPieceTraj.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkEnv/TrkEnv.hh"
#include "TrkBase/HelixTraj.hh"
#include "CLHEP/Geometry/HepPoint.h"
#include "difAlgebra/DifPoint.hh"
#include "difAlgebra/DifVector.hh"
#include "difAlgebra/DifIndepPar.hh"

//------------------------------------------------------------------------
KalPairConduit::KalPairConduit(const TrkDifPieceTraj& plusTraj,
			       const TrkDifPieceTraj& minusTraj,
			       const BbrVectorErr& beamMom) :
  _nreps(0), _nsites(0),  _niter(0),
  _plusFltLen(), _minusFltLen(),
  _prodPoint(), _prodPointChi2(0.), _converged(false), _needsUpdate(true),
  _beamMom(beamMom)
{
  _reps[plusRep] = _reps[minusRep] = 0;
  _sites[plusRep] = _sites[minusRep] = 0;
  // Set the LookAtMe state for each Site to true
  _siteLAM[plusRep] = _siteLAM[minusRep] = true;
  _constraintPar[plusRep] = new KalParams(5);
  _constraintPar[minusRep] = new KalParams(5);

  TrkPoca poca(plusTraj, 0.0, minusTraj, 0.0);
  _fltlens[plusRep] = poca.flt1();
  _fltlens[minusRep] = poca.flt2();

  double dummy;
  TrkParams tempPar =
    *plusTraj.localTrajectory(poca.flt1(), dummy)->parameters();
  _inwardPar[plusRep] = new KalParams(tempPar);
  tempPar = *minusTraj.localTrajectory(poca.flt2(), dummy)->parameters();
  _inwardPar[minusRep] = new KalParams(tempPar);

  updateConstraints();
}
//------------------------------------------------------------------------
KalPairConduit::~KalPairConduit()
{
  delete _constraintPar[plusRep];
  delete _constraintPar[minusRep];
  delete _inwardPar[plusRep];
  delete _inwardPar[minusRep];
}
//------------------------------------------------------------------------
TrkErrCode KalPairConduit::coordinateFit()
{
  // The procedure is to perform the individual Rep fits using the existing
  // constraint information.  Then, update the constraint information
  // using the inward results from the fits.  Iterate this until the
  // constraints stop changing (or we've run out of iterations).

  // Should check here if the inward fit is valid.  If it is, no need to
  // change the constraint site.  If it isn't, do the inward fit, update
  // the constraint site, and continue.  Inward fit could be invalid if
  // someone disables a HoT, for example.
  // This won't produce a completely valid fit, but will invoke the
  // pairSite's process method.  We'll tell it to upload its constraint info
  TrkErrCode fitState;
  if (_reps[plusRep]->needsFit(trkIn)) {
    fitState = _reps[plusRep]->KalRep::fit(trkIn);
    if (fitState.failure()) return fitState;
  }
  if (_reps[minusRep]->needsFit(trkIn)) {
    fitState = _reps[minusRep]->KalRep::fit(trkIn);
    if (fitState.failure()) return fitState;
  }
  // If either of the inward fit results have changed, compute the constraint
  // information.  This needs to be done for both reps even if only one has
  // changed. The _needsUpdate flag will be set via the PairSite's process
  // method if the new inward fit is different than the currently cached one.
  if (_needsUpdate) updateConstraints();

  // Now we are ready to do a full blown fit.  The PairSites will retrieve
  // their constraints from this conduit when their update method is invoked.
  fitState = _reps[plusRep]->KalRep::fit();
  _reps[plusRep]->addHistory(fitState, "KalPairFit");
  if (fitState.failure()) return fitState;
  // Make it fail if it didn't converge
//   if (2 == fitState.success()) {
//     fitState.setFailure(2);
//     return fitState;
//   }

  fitState = _reps[minusRep]->KalRep::fit();
  _reps[minusRep]->addHistory(fitState, "KalPairFit");
  if (fitState.failure()) return fitState;
//   if (2 == fitState.success()) {
//     fitState.setFailure(2);
//     return fitState;
//   }

  calcProdPoint();

  return fitState;
}
//------------------------------------------------------------------------
KalParams KalPairConduit::getConstraint(KalPairRep* rep)
{
  int index = thisRepIndex(rep);
  _siteLAM[index] = false;
  return *(_constraintPar[index]);
}
//------------------------------------------------------------------------
KalParams KalPairConduit::getConstraint(KalPairSite* site)
{
  return getConstraint(_reps[thisSiteIndex(site)]);
}
//------------------------------------------------------------------------
double KalPairConduit::getFltLen(KalPairRep* rep)
{
  int index = thisRepIndex(rep);
  return _fltlens[index];
}
//------------------------------------------------------------------------
double KalPairConduit::getFltLen(KalPairSite* site)
{
  return getFltLen(_reps[thisSiteIndex(site)]);
}
//------------------------------------------------------------------------
void KalPairConduit::killedBy(KalPairRep* thisRep)
{
  // Tell the other rep to set its pointer to me to NULL
  _reps[otherRepIndex(thisRep)]->nullConduit();
}
//------------------------------------------------------------------------
KalPairRep* KalPairConduit::otherPairRep(KalPairRep* rep)
{
  return _reps[otherRepIndex(rep)];
}
//------------------------------------------------------------------------
int KalPairConduit::otherRepIndex(KalPairRep* thisRep)
{
  // Return the index of the other KalPairRep
  return (1 - thisRepIndex(thisRep));
}
//------------------------------------------------------------------------
int KalPairConduit::thisRepIndex(KalPairRep* thisRep)
{
  // Find the index of this rep
  for (int i = 0; i < 2; i++)
    if (_reps[i] == thisRep)
      return i;

  // If we get here, something is horribly wrong!
  assert(true);
  return -1;
}
//------------------------------------------------------------------------
int KalPairConduit::thisSiteIndex(KalPairSite* thisSite)
{
  // Find the index of this site
  for (int i = 0; i < 2; i++)
    if (_sites[i] == thisSite)
      return i;

  // If we get here, something is horribly wrong!
  assert(true);
  return -1;
}
//------------------------------------------------------------------------
void KalPairConduit::addRep(KalPairRep* newRep)
{
  // Add this rep to the rep list
  assert(_nreps < 2);
  int index;
  if (newRep->charge() > 0)
    index = plusRep;
  else
    index = minusRep;

  assert(0 == _reps[index]);
  _reps[index] = newRep;
  _nreps++;

  return;
}
//------------------------------------------------------------------------
void KalPairConduit::addSite(KalPairSite* newSite, KalPairRep* rep)
{
  // Add this site with the same index as the rep (just to avoid confusion)
  int index = thisRepIndex(rep);
  // Being pretty assertive here...
  assert(index < 2);
  assert(_sites[index] == 0);
  assert(_nsites < 2);
  _sites[index] = newSite;
  _nsites++;
  return;
}
//------------------------------------------------------------------------
bool KalPairConduit::checkLAM(KalPairSite* site)
{
  return _siteLAM[thisSiteIndex(site)];
}
//------------------------------------------------------------------------
void KalPairConduit::uploadParams(const KalParams& params, KalPairSite* site)
{
  // Check if these parameters that the site is trying to upload are
  // 'significantly' different than the ones already here.  If they are,
  // replace the exisiting parameters and set a flag to let this guy know
  // that it should update the constraints.  Shouldn't perform the update
  // yet, because the other site may still have to upload its parameters.
  // If existing params aren't valid, then just do the update.

  double chisq = 9999.;
  if (_inwardPar[thisSiteIndex(site)]->matrixOK())
    chisq = params.chisq(*_inwardPar[thisSiteIndex(site)]);
//  cout << "In KalPairConduit: chisq = " << chisq << endl;
  if (chisq > 0.000001) {
    _needsUpdate = true;
    *_inwardPar[thisSiteIndex(site)] = params;
  }
  
  return;
}
//------------------------------------------------------------------------
bool KalPairConduit::updateConstraints()
{
  // This should take the untransformed parameters and construct constraints
  // to be used by the other rep.  Should also tell the pairSites that
  // they have new constraints to download.

  // Make HelixTrajs out of the untransformed parameters...
  TrkExchangePar plusPar(_inwardPar[plusRep]->parameterVector(),
			 _inwardPar[plusRep]->covarianceMatrix());
  TrkExchangePar minusPar(_inwardPar[minusRep]->parameterVector(),
			  _inwardPar[minusRep]->covarianceMatrix());
  HelixTraj plusTraj(plusPar);
  HelixTraj minusTraj(minusPar);

  KalPairVisitor plusVis(plusTraj);
  KalPairVisitor minusVis(minusTraj);

  // Determine the flightlength at which to perform the transformation
  _plusFltLen =
    plusVis.fitFlightLength(0., _beamMom,
			    _inwardPar[minusRep]->trackParameters());
  _minusFltLen =
    minusVis.fitFlightLength(0., _beamMom,
			     _inwardPar[plusRep]->trackParameters());

  // Do the transformation
  TrkParams plusConPar(5), minusConPar(5);
  bool isOk = (plusVis.transformParams(_plusFltLen, _beamMom, minusConPar) &&
	       minusVis.transformParams(_minusFltLen, _beamMom, plusConPar));
  if (!isOk) return false;

  *_constraintPar[plusRep] = plusConPar;
  _siteLAM[plusRep] = true;
  // This method gets called in the ctor, before the _sites exist.
  // Do I even need to do this?
  // Now actually call the update method of the KalPairSites
  if (0 != _sites[plusRep]) {
    _sites[plusRep]->invalidateSite(trkIn);
    _sites[plusRep]->invalidateSite(trkOut);
    if ( 0 != _reps[plusRep] ) {
      _reps[plusRep]->setCurrent(false);
      _reps[plusRep]->resetFit();
      _sites[plusRep]->update(_reps[plusRep]->referenceTraj(),
			      _reps[plusRep]->refMomentum());
    } 
  }

  *_constraintPar[minusRep] = minusConPar;
  _siteLAM[minusRep] = true;
  if (0 != _sites[minusRep]) {
    _sites[minusRep]->invalidateSite(trkIn);
    _sites[minusRep]->invalidateSite(trkOut);
    if ( 0 != _reps[minusRep] ) {
      _reps[minusRep]->setCurrent(false);
      _reps[minusRep]->resetFit();
      _sites[minusRep]->update(_reps[minusRep]->referenceTraj(),
			      _reps[minusRep]->refMomentum());
    } 
  }

  _needsUpdate = false;

  return true;
}
//------------------------------------------------------------------------
void KalPairConduit::calcProdPoint()
{
  // The flight lengths calculated by the visitors can be used to determine
  // the production point that minimizes the chi2.  This is potentially
  // different than the point determined by the poca.  By taking a weighted
  // average of the space point defined by the inwards fits at this 
  // flight length, without propagating errors on the flight length, we
  // should get two uncorrelated measures of the prodPoint (I think).

  // First, make a HelixTraj out of the inward fit parameters
  HelixTraj plusTraj(_inwardPar[plusRep]->parameterVector(),
		     _inwardPar[plusRep]->covarianceMatrix());
  HelixTraj minusTraj(_inwardPar[minusRep]->parameterVector(),
		     _inwardPar[minusRep]->covarianceMatrix());

  DifPoint posD;
  DifVector dirD;
  plusTraj.getDFInfo2(_plusFltLen.value(), posD, dirD);
  HepMatrix errPlus = posD.errorMatrix(posD.x.indepPar()->covariance());
  HepPoint pointPlus(posD.x.number(), posD.y.number(), posD.z.number());
  BbrError symErrPlus(3);
  symErrPlus.assign(errPlus);
  const BbrPointErr xPlus(pointPlus, symErrPlus);
  
  minusTraj.getDFInfo2(_minusFltLen.value(), posD, dirD);
  HepMatrix errMinus = posD.errorMatrix(posD.x.indepPar()->covariance());
  HepPoint pointMinus(posD.x.number(), posD.y.number(), posD.z.number());
  BbrError symErrMinus(3);
  symErrMinus.assign(errMinus);
  const BbrPointErr xMinus(pointMinus, symErrMinus);

  //  const BbrPointErr xPlus(_reps[plusRep]->positionErr(_plusFltLen.value()));
  //  const BbrPointErr xMinus(_reps[minusRep]->positionErr(_minusFltLen.value()));
  HepSymMatrix v1(xPlus.covMatrix());
  HepSymMatrix v2(xMinus.covMatrix());
  int ierr;
  v1.invert(ierr);
  if (0 != ierr) return;
  v2.invert(ierr);
  if (0 != ierr) return;
  HepSymMatrix vsum = v1 + v2;
  vsum.invert(ierr);
  if (0 != ierr) return;
  HepVector x1(3);
  x1[0] = xPlus.x(); x1[1] = xPlus.y(); x1[2] = xPlus.z();
  HepVector x2(3);
  x2[0] = xMinus.x(); x2[1] = xMinus.y(); x2[2] = xMinus.z();
  HepVector xAvg = vsum*(v1*x1 + v2*x2);
  HepPoint x(xAvg[0], xAvg[1], xAvg[2]);
  BbrPointErr newProd(x, vsum);
  _prodPoint = newProd;
  _prodPointChi2 = v1.similarity(x1-xAvg) + v2.similarity(x2-xAvg);
  
  return;
}
//------------------------------------------------------------------------
