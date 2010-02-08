//--------------------------------------------------------------------------
// File and Version Information:
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
#include "BbrGeom/BbrPointErr.hh"
#include "BbrGeom/BbrVectorErr.hh"
#include "BField/BField.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "KalmanTrack/KalPairVisitor.hh"
#include "TrkEnv/TrkEnv.hh"
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "TrkBase/TrkParams.hh"
#include "TrkBase/TrkPoca.hh"

//------------------------------------------------------------------------
KalPairVisitor::~KalPairVisitor() {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
KalPairVisitor::KalPairVisitor(const TrkSimpTraj& theTraj) {
//------------------------------------------------------------------------

  theTraj.visitAccept(this);
}

//------------------------------------------------------------------------
void
KalPairVisitor::trkVisitHelixTraj(const HelixTraj* theTraj) {
//------------------------------------------------------------------------
// set the "array"

  _ht = theTraj;
  _ct = 0;
  _nt = 0;
  _lt = 0;
}

//------------------------------------------------------------------------
void
KalPairVisitor::trkVisitCircleTraj(const TrkCircleTraj* theTraj) {
//------------------------------------------------------------------------
// set the "array"

  _ht = 0;
  _ct = theTraj;
  _nt = 0;
  _lt = 0;
}

//------------------------------------------------------------------------
void
KalPairVisitor::trkVisitNeutTraj(const NeutTraj* theTraj) {
//------------------------------------------------------------------------
// set the "array"

  _ht = 0;
  _ct = 0;
  _nt = theTraj;
  _lt = 0;
}

//------------------------------------------------------------------------
void
KalPairVisitor::trkVisitLineTraj(const TrkDifLineTraj* theTraj) {
//------------------------------------------------------------------------
// set the "array"

  _ht = 0;
  _ct = 0;
  _nt = 0;
  _lt = theTraj;
}
//------------------------------------------------------------------------
BbrDoubleErr
KalPairVisitor::fitFlightLength(double fltlen, BbrVectorErr mom,
				TrkParams otherParams,
				int level) const {

  assert (0 != _ht);

  level++;
  // This duplicates a lot of code in transformParams.  Should probably
  // do some redesign at some point.

  // Use this helix to determine the momentum and position
  const BField* field = gblEnv->getTrk()->magneticField();
  Hep3Vector oldMom =
    TrkMomCalculator::vecMom(*_ht, *field, fltlen);
  int charge = TrkMomCalculator::charge(*_ht, *field, fltlen);
  HepPoint pos = _ht->position(fltlen);

  // The new momentum is the constrained mom minus the old mom
  Hep3Vector newMom = (Hep3Vector)mom - oldMom;

  // Now, get new parameters based on the new momentum
  int newCharge = -1*charge;
  TrkExchangePar newTrkPar =
    TrkHelixUtils::helixFromMom(pos, newMom, newCharge, *field);

  // Transform covariance matrix
  // First, construct the Jacobian for the parameter transformation.
  // For notation clarification, 'eta' refers to vector of parameters,
  // 'q' refers to momentum and position, {p, x}. Transformed (new) quantites
  // are primed.
  // Jacobian is d(eta')/d(eta) = [d(eta')/d(q')][d(q')/d(eta)]
  
  // Get matrix of derivatives of new parameters wrt momentum and position.
  HepMatrix DetaDq = derivsParWrtCoord(newMom, pos, newCharge);
  HepMatrix DqDeta = derivsCoordWrtPar(_ht->parameters()->parameter(),
				       charge, fltlen);
  // Jacobian is Deta'Dq' times Dq'Deta
  HepMatrix JacEta = DetaDq * DqDeta;

  // New covariance matrix is V' = JacEta * V * JacEta.T
  HepSymMatrix newCov = _ht->parameters()->covariance();
  newCov = newCov.similarity(JacEta);

  // Now added error on momentum constraint to covariance matrix
  HepMatrix JacBeam = DetaDq * derivsCoordWrtBeam();
  newCov += mom.covMatrix().similarity(JacBeam);

  // Now determine the difference between the transformed parameters
  // and the passed parameters.  We want to find the flightlength that
  // minimizes this difference.

  HepVector deltaPar = otherParams.parameter() -
    newTrkPar.params();
  HepSymMatrix deltaCov = otherParams.covariance() + newCov;
    //    newTrkPar.covariance();
  int ifail;
  HepSymMatrix deltaCovInv = deltaCov.inverse(ifail);

  // Now we just need the derivatives of the transformed parameters
  // w.r.t. the flightlength

  // For notation clarification, 'eta' refers to vector of parameters,
  // 'q' refers to momentum and position, {p, x}. Transformed (new) quantites
  // are primed.
  // Jacobian is d(eta')/d(eta) = [d(eta')/d(q')][d(q')/d(eta)]
  
  // Get matrix of derivatives of new parameters wrt momentum and position.
  //  HepMatrix DetaDq = derivsParWrtCoord(newMom, pos, newCharge);
  HepVector JacFlt = DetaDq*derivsCoordWrtFlt(_ht->parameters()->parameter(),
					      charge, fltlen);

  // This should just be a 1x1 matrix
  double sig2Inv = deltaCovInv.similarity(JacFlt);
  double sig2 = 1./sig2Inv;
  HepVector deltaFltLen = sig2*JacFlt.T()*deltaCovInv*deltaPar;
  double delta = deltaFltLen[0];
  //  chi2 = deltaCovInv.similarity(deltaPar);
  //  cout << " Level in fitFlightLength = " << level << endl;
  //  cout << " Chi2 in fitFlightLength = " << chi2 << endl;

  // Should probably iterate this fit, since the transformed covariance
  // matrix also depends on the flightlength.
  if( fabs(delta) < 0.00001 || level > 5) {
    BbrDoubleErr returnVal(fltlen+delta, sig2);
    return returnVal;
  } else return fitFlightLength(fltlen+delta, mom, otherParams, level);
}
//------------------------------------------------------------------------
bool
KalPairVisitor::transformParams(BbrDoubleErr fltlen, BbrVectorErr mom,
				TrkParams& params) const {
//------------------------------------------------------------------------
  // This only knows how to handle helix trajectories
  if (0 == _ht) return false;

  // Use this helix to determine the momentum and position
  const BField* field = gblEnv->getTrk()->magneticField();
  Hep3Vector oldMom =
    TrkMomCalculator::vecMom(*_ht, *field, fltlen.value());
  int charge = TrkMomCalculator::charge(*_ht, *field, fltlen.value());
  HepPoint pos = _ht->position(fltlen.value());

  // The new momentum is the constrained mom minus the old mom
  Hep3Vector newMom = (Hep3Vector)mom - oldMom;

  // Now, get new parameters based on the new momentum
  int newCharge = -1*charge;
  TrkExchangePar newTrkPar =
    TrkHelixUtils::helixFromMom(pos, newMom, newCharge, *field);

  // Transform covariance matrix
  // First, construct the Jacobian for the parameter transformation.
  // For notation clarification, 'eta' refers to vector of parameters,
  // 'q' refers to momentum and position, {p, x}. Transformed (new) quantites
  // are primed.
  // Jacobian is d(eta')/d(eta) = [d(eta')/d(q')][d(q')/d(eta)]
  
  // Get matrix of derivatives of new parameters wrt momentum and position.
  HepMatrix DetaDq = derivsParWrtCoord(newMom, pos, newCharge);
  HepMatrix DqDeta = derivsCoordWrtPar(_ht->parameters()->parameter(),
				       charge, fltlen.value());
  // Jacobian is Deta'Dq' times Dq'Deta
  HepMatrix JacEta = DetaDq * DqDeta;

  // New covariance matrix is V' = JacEta * V * JacEta.T
  HepSymMatrix newCov = _ht->parameters()->covariance();
  newCov = newCov.similarity(JacEta);

  // Now added error on momentum constraint to covariance matrix
  HepMatrix JacBeam = DetaDq * derivsCoordWrtBeam();
  newCov += mom.covMatrix().similarity(JacBeam);

  // Add error on flight length position of pair site
  HepVector JacFlt = DetaDq*derivsCoordWrtFlt(_ht->parameters()->parameter(),
 					      charge, fltlen.value());
  newCov += fltlen.covariance()*vT_times_v(JacFlt);

  // Construct new TrkParams
  TrkParams newParams(newTrkPar.params(), newCov);

  params = newParams;

  return true;
}
//------------------------------------------------------------------------
const HepMatrix
KalPairVisitor::derivsCoordWrtPar(TrkExchangePar params, int charge,
				  double fltlen) const {
//------------------------------------------------------------------------
  assert(0 != _ht);  // Only works for helices

  HepMatrix derivs(nCoord, TrkExchangePar::nParam, 0);

  // Calcualte some common stuff
  const BField* field = gblEnv->getTrk()->magneticField();
  const double A = -BField::cmTeslaToGeVc*field->bFieldNominal()*charge;
  const double d0 = params.d0();
  const double phi0 = params.phi0();
  const double omega = params.omega();
  const double tanDip = params.tanDip();
  const double cosDip = 1./sqrt(1. + tanDip*tanDip);
  const double sinDip = cosDip*tanDip;
  const double translen = fltlen*cosDip;
  const double arc = translen*omega;
  const double cosPhi = cos(phi0);
  const double sinPhi = sin(phi0);
  const double cosAng = cos(phi0 + arc);
  const double sinAng = sin(phi0 + arc);
  const double omega2 = sqr(omega);

  assert (0 != omega);

  // pX wrt phi0
  derivs[pX][TrkExchangePar::ex_phi0] = (A/omega)*sinAng;
  // pX wrt omega
  derivs[pX][TrkExchangePar::ex_omega] = (A/omega2)*(cosAng + arc*sinAng);
  // pX wrt tanDip
  derivs[pX][TrkExchangePar::ex_tanDip] = -A*sinAng*translen*sinDip*cosDip;

  // pY wrt phi0
  derivs[pY][TrkExchangePar::ex_phi0] = -(A/omega)*cosAng;
  // pY wrt omega
  derivs[pY][TrkExchangePar::ex_omega] = (A/omega2)*(sinAng - arc*cosAng);
  // pY wrt tanDip
  derivs[pY][TrkExchangePar::ex_tanDip] = A*cosAng*translen*sinDip*cosDip;

  // pZ wrt omega
  derivs[pZ][TrkExchangePar::ex_omega] = (A/omega2)*tanDip;
  // pZ wrt tanDip
  derivs[pZ][TrkExchangePar::ex_tanDip] = -A/omega;

  // X wrt d0
  derivs[X][TrkExchangePar::ex_d0] = -sinPhi;
  // X wrt phi0
  derivs[X][TrkExchangePar::ex_phi0] = -d0*cosPhi + (cosAng - cosPhi)/omega;
  // X wrt omega
  derivs[X][TrkExchangePar::ex_omega] = (1./omega2) * 
    (sinPhi - sinAng + arc*cosAng);
  // X wrt tanDip
  derivs[X][TrkExchangePar::ex_tanDip] = -cosAng*translen*sinDip*cosDip;

  // Y wrt d0
  derivs[Y][TrkExchangePar::ex_d0] = cosPhi;
  // Y wrt phi0
  derivs[Y][TrkExchangePar::ex_phi0] = -d0*sinPhi + (sinAng - sinPhi)/omega;
  // Y wrt omega
  derivs[Y][TrkExchangePar::ex_omega] = (1./omega2) * 
    (-cosPhi + cosAng + arc*sinAng);
  // Y wrt tanDip
  derivs[Y][TrkExchangePar::ex_tanDip] = -sinAng*translen*sinDip*cosDip;

  // Z wrt z0
  derivs[Z][TrkExchangePar::ex_z0] = 1.;
  // Z wrt tanDip
  derivs[Z][TrkExchangePar::ex_tanDip] = translen*cosDip*cosDip;

  return derivs;
}
//------------------------------------------------------------------------
const HepMatrix
KalPairVisitor::derivsParWrtCoord(Hep3Vector mom, HepPoint pos,
				  int charge) const {
//------------------------------------------------------------------------
  assert (0 != _ht);  // Only works for helices

  HepMatrix derivs(TrkExchangePar::nParam, nCoord, 0);

  // Calcualte some common stuff
  const BField* field = gblEnv->getTrk()->magneticField();
  const double A = -BField::cmTeslaToGeVc*field->bFieldNominal()*charge;
  const double pT = mom.perp();
  const double pT2 = mom.perp2();
  const double pT3 = pow(pT, 3.);
  const double pYx = mom.y() - A*pos.x();
  const double pXy = mom.x() + A*pos.y();
  const double denom2 = sqr(pYx) + sqr(pXy);
  const double denom = sqrt(denom2);
  const double r2 = sqr(pos.x()) + sqr(pos.y());

  assert(0. != A);
  assert(0. != pT);
  assert(0. != denom);

  // d0 wrt pX
  derivs[TrkExchangePar::ex_d0][pX] = (1/A) * (pXy/denom - mom.x()/pT);
  // d0 wrt pY
  derivs[TrkExchangePar::ex_d0][pY] = (1/A) * (pYx/denom - mom.y()/pT);
  // d0 wrt X
  derivs[TrkExchangePar::ex_d0][X] = -pYx/denom;
  // d0 wrt Y
  derivs[TrkExchangePar::ex_d0][Y] = pXy/denom;

  // phi0 wrt pX
  derivs[TrkExchangePar::ex_phi0][pX] = -pYx/denom2;
  // phi0 wrt pY
  derivs[TrkExchangePar::ex_phi0][pY] = pXy/denom2;
  // phi0 wrt X
  derivs[TrkExchangePar::ex_phi0][X] = -A * pXy/denom2;
  // phi0 wrt Y
  derivs[TrkExchangePar::ex_phi0][Y] = -A * pYx/denom2;

  // omega wrt pX
  derivs[TrkExchangePar::ex_omega][pX] = -A * mom.x()/pT3;
  // omega wrt pY
  derivs[TrkExchangePar::ex_omega][pY] = -A * mom.y()/pT3;

  // z0 wrt pX
  derivs[TrkExchangePar::ex_z0][pX] = 
    mom.z()*(sqr(mom.x())*pos.x() + 2*mom.x()*mom.y()*pos.y() -
	     sqr(mom.y())*pos.x() + A*mom.y()*r2) / (pT2 * denom2);
  // z0 wrt pY
  derivs[TrkExchangePar::ex_z0][pY] = 
    -mom.z()*(sqr(mom.x())*pos.y() - 2*mom.x()*mom.y()*pos.x() -
	      sqr(mom.y())*pos.y() + A*mom.x()*r2) / (pT2 * denom2);
  // z0 wrt pZ
  derivs[TrkExchangePar::ex_z0][pZ] =
    -(1./A) * atan(A*(mom.x()*pos.x() + mom.y()*pos.y())/
		   (mom.x()*pXy + mom.y()*pYx));
  // z0 wrt X
  derivs[TrkExchangePar::ex_z0][X] = -mom.z()*pXy/denom2;
  // z0 wrt Y
  derivs[TrkExchangePar::ex_z0][Y] = -mom.z()*pYx/denom2;
  // z0 wrt Z
  derivs[TrkExchangePar::ex_z0][Z] = 1.;

  // tanDip wrt pX
  derivs[TrkExchangePar::ex_tanDip][pX] = -mom.x()*mom.z()/pT3;
  // tanDip wrt pY
  derivs[TrkExchangePar::ex_tanDip][pY] = -mom.y()*mom.z()/pT3;
  // tanDip wrt pZ
  derivs[TrkExchangePar::ex_tanDip][pZ] = 1./pT;


  return derivs;
}
//------------------------------------------------------------------------
const HepMatrix
KalPairVisitor::derivsCoordWrtBeam() const {
//------------------------------------------------------------------------
  assert (0 != _ht);  // Only works for helices

  HepMatrix derivs(nCoord, Hep3Vector::NUM_COORDINATES, 0);

  // Very simple: q' = q + p_Beam.  We can ignore sign of p_Beam (it is
  // different for the inward and outward fit) because it will cancel during
  // the similarity transformation

  derivs[pX][Hep3Vector::X] = 1.;
  derivs[pY][Hep3Vector::Y] = 1.;
  derivs[pZ][Hep3Vector::Z] = 1.;

  return derivs;
}
//------------------------------------------------------------------------
const HepVector
KalPairVisitor::derivsCoordWrtFlt(TrkExchangePar params, int charge,
				  double fltlen) const {
//------------------------------------------------------------------------
  assert (0 != _ht);  // Only works for helices

  HepVector derivs(nCoord, 0);

  // Calcualte some common stuff
  const BField* field = gblEnv->getTrk()->magneticField();
  const double A = -BField::cmTeslaToGeVc*field->bFieldNominal()*charge;
  const double phi0 = params.phi0();
  const double tanDip = params.tanDip();
  const double cosDip = 1./sqrt(1. + sqr(tanDip));
  const double arc = fltlen*cosDip*params.omega();
  const double cosAng = cos(phi0 + arc);
  const double sinAng = sin(phi0 + arc);

  // pX wrt translen
  derivs[pX] = A*sinAng;
  // pY wrt translen
  derivs[pY] = -A*cosAng;

  // X wrt translen
  derivs[X] = cosAng;
  // Y wrt translen
  derivs[Y] = sinAng;
  // Z wrt translen
  derivs[Z] = tanDip;

  // dervis wrt fltlen
  derivs *= cosDip;

  return derivs;
}
//------------------------------------------------------------------------
double
KalPairVisitor::calcSigma2FltLen(const double fltIn,
				 const BbrPointErr x) const
{
  assert (0 != _ht);

  // This does nasty numeric derivative calculations making a billion
  // calls to poca.

  // First, determine step sizes to make for parameter variation.  Will be
  // a constant fraction of the parameter error.
  const double sigFrac(0.1);
  // Get matrix of derivatives of fltlen w.r.t. direct parameter variation
  HepMatrix dfdpar(1, HelixTraj::NHLXPRM, 0);
  int i;
  for (i = 0; i < HelixTraj::NHLXPRM; i++) {
    double parStep = sigFrac*sqrt(_ht->parameters()->covariance()[i][i]);
    // Make a new helix trajectory with the tweaked parameter
    TrkExchangePar newpar(_ht->parameters()->parameter(),
			  _ht->parameters()->covariance());
    newpar.params()[i] += parStep;
    HelixTraj newTraj(newpar, -99999., 99999., _ht->referencePoint());
    TrkPoca poca(newTraj, fltIn, x);
    assert (0. != parStep);
    dfdpar[0][i] = (poca.flt1() - fltIn)/parStep;
  }

  // Get matrix of derivatives w.r.t. point error
  HepMatrix dfdx(1, 3, 0);
  for (i = 0; i < 3; i++) {
    double xStep = sigFrac*sqrt(x.covMatrix()[i][i]);
    //    cout << " i, xStep = " << i << "  " << xStep << endl;
    Hep3Vector delX;
    delX[i] = xStep;
    assert(0. != xStep);
    HepPoint newX(x);
    newX += delX;
    TrkPoca poca(*_ht, fltIn, newX);
    dfdx[0][i] = (poca.flt1() - fltIn)/xStep;
  }

  HepSymMatrix temp = _ht->parameters()->covariance().similarity(dfdpar);
  temp += x.covMatrix().similarity(dfdx);

  return temp[0][0];
}
//------------------------------------------------------------------------




