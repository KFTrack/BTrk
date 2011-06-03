//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHelixUtils.cc,v 1.6 2008/04/28 06:01:53 brownd Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner, code stolen from TrkParms class
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/NeutParams.hh"
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <math.h>
#include "TrkBase/TrkExchangePar.hh"
#include "BField/BField.hh"
#include "BbrGeom/BbrAngle.hh"
#include "BbrGeom/BbrPointErr.hh"
#include "BbrGeom/BbrVectorErr.hh"
#include "difAlgebra/DifNumber.hh"
#include "difAlgebra/DifArray.hh"
#include "BaBar/Constants.hh"

const double TrkHelixUtils::small(1e-14);

//------------------------------------------------------------------------
HepMatrix TrkHelixUtils::jacobianExtrapolate(const TrkExchangePar& par, 
					  double fltNew) {
//------------------------------------------------------------------------
  
//----------------------------------------------------------------------
//      Compute and return the jacobian that takes the covariance matrix
//	from fltOld to fltNew
//
//      "fltLen" -- the signed 3-d pathlength a particle travels 
//                     along the orbit starting from the point on the orbit 
//                     close to the origin.
//
//      Each element in this matrix is a partial derivative of an orbit
//      parameter at some reference point to an orbit parameter at
//      the fit point.  What is kept fixed in taking these partial 
//      derivatives is that the orbit parameters are those at the point
//      on the orbit that is closest in the x-y plane to the reference point.
//      This transform matrix has the property that transform(-arclength) 
//      is the inverse of transform(arclength).
//      (repaired by Gerry Lynch, I think -- sfs)
//----------------------------------------------------------------------
//      A more general calculation of this tranformation matrix is in
//      TrkBase/HelixTraj::paramFunc.  jacobianExtrapolate assumes that
//      the reference points are on the track.  paramFunc does not.
//----------------------------------------------------------------------

  HepMatrix transform(5, 5, 0);

  double tanDip = par.tanDip();
  double omega = par.omega();
  // Convert to 2-d arclength
  double darc = fltNew / sqrt(1. + tanDip * tanDip);
  double dphi    = omega * darc;
  double cosDphi = cos(dphi);
  double sinDphi = sin(dphi);

  // partials of d0
 
  transform[0][0] = cosDphi;
  transform[0][1] = sinDphi / omega;
  transform[0][2] = (1.0-cosDphi) / (omega * omega);

  // partials of phi0  

  transform[1][0] = -sinDphi * omega;
  transform[1][1] =  cosDphi;
  transform[1][2] =  sinDphi / omega;

  // partials of omega  

  transform[2][2] = 1.0;

  // partials of z0 
 
  transform[3][0] = -tanDip * sinDphi;
  transform[3][1] = -tanDip * (1.0-cosDphi) / omega;
  transform[3][2] = -tanDip * (dphi-sinDphi) / (omega*omega);
  transform[3][3] =  1.;
  transform[3][4] =  darc;

  // partials of tanDip

  transform[4][4] = 1.;

  return transform;
}


//----------------------------------------------------------------------
HepSymMatrix TrkHelixUtils::extrapolateCov(TrkExchangePar& par, 
					   double fltNew) {
//----------------------------------------------------------------------

  return par.covariance().similarity(jacobianExtrapolate(par, fltNew));  
}

TrkExchangePar TrkHelixUtils::helixFromMom(const HepPoint& pos, 
       	      const Hep3Vector& pmom, double sign, const BField& fieldmap) {
	static HepVector pars(5);
	static double fltlen;
	helixFromMom(pars,fltlen,pos,pmom,sign,fieldmap);
	return TrkExchangePar(pars);
}

//----------------------------------------------------------------------
void TrkHelixUtils::helixFromMom(HepVector& pars, double& fltlen, const HepPoint& pos, 
const Hep3Vector& pmom, double sign, const BField& fieldmap) {
  helixFromMom(pars,fltlen,pos,pmom,sign,fieldmap.bFieldNominal());
}


//----------------------------------------------------------------------
void TrkHelixUtils::helixFromMom(HepVector& pars, double& fltlen, const HepPoint& pos, 
const Hep3Vector& pmom, double sign, double Bval) {
//----------------------------------------------------------------------
  // Before September 2005 this used the equations
  //As documented in 
  //R.~Luchsinger and C.~Grab, Comp.~Phys.~Comm.~{\bf 76}, 263-280 (1993).
  //Equation 14 on page 270.
  //Note: We use the opposite sign for d0.
  //We use tandip and not the dip angle itself.
  // The code that was in this routine before September got the wrong
  //  answer when the turning angle from the beginning of the track
  //  and the correct reference point is greater than 90 degree.  

  register double phip,rho,pt;
  
  double px = pmom.x();
  double py = pmom.y();
  pt=sqrt(px*px+py*py);
  if (pt < small)   pt = small;  // hack to avoid pt=0 tracks
  if (fabs(px) < small) px = (px<0.0) ? -small : small; // hack to avoid pt=0 tracks
    
  pars(3)=-BField::cmTeslaToGeVc*Bval*sign/pt;  //omega
  pars(5)=pmom.z()/pt;  //tandip
  rho=1./pars(3);
//  double radius = fabs(rho);
  phip=atan2(py,px);
//  double cosphip=px/pt; double sinphip=py/pt; // this should be more efficient; test!
  double cosphip=cos(phip);
  double sinphip=sin(phip);
  double cx = pos.x()-rho*sinphip;
  double cy = pos.y()+rho*cosphip;
  double phi0 = atan2(-cx,cy);
  double sphi0 = sin(phi0);
  double cphi0 = cos(phi0);
  if( (cy*cphi0-cx*sphi0)*rho < 0.0 ){
  // wrong angular momentum: fix
    phi0 = BbrAngle(phi0+Constants::pi).rad();
    sphi0 = -sphi0;
    cphi0 = -cphi0;
  }
//  pars(2) = phi0;
  double d0;
  if(fabs(sphi0)>0.5)
    d0 = -cx/sphi0-rho;
  else
    d0 = cy/cphi0-rho;
  pars(1) = d0;
  pars(2) = phi0;
  double translen = rho*BbrAngle(phip-pars(2)).rad();
  pars(4) = pos.z()-translen*pars(5); // z0
  double oldflt = fltlen;
// flightlength; measured in 3-space!!!
  double ffact = sqrt(1.0+pars(5)*pars(5));
  fltlen = translen*ffact;
// adjust flightlength and z0 to correspond to the same # of loops as the input
  if(fabs(fltlen - oldflt) > fabs(rho)){
    double floop = rho*Constants::twoPi*ffact;
    int nloop = (int)rint((oldflt-fltlen)/floop);
    fltlen += nloop*floop;
    pars(4) -= nloop*rho*Constants::twoPi*pars(5);
  }
}

//----------------------------------------------------------------------
TrkExchangePar TrkHelixUtils::helixFromMomErr(const BbrPointErr& pos,
              const BbrVectorErr& pmom,const HepMatrix& cxp, double sign, 
              const BField& fieldmap) {
//----------------------------------------------------------------------

  DifNumber xDF(pos.x(),1,6), yDF(pos.y(),2,6), zDF(pos.z(),3,6);
  DifNumber pxDF(pmom.x(), 4, 6);
  DifNumber pyDF(pmom.y(), 5, 6);
  DifNumber pzDF(pmom.z(), 6, 6);

  static DifNumber phip, cosphip, sinphip, gamma;
  static DifArray pars(5,6);

  DifNumber invpt = pxDF;
  invpt *= pxDF;
  invpt += pyDF*pyDF;
  invpt = sqrt(invpt);

  if (invpt < small) invpt = small;  // hack to avoid pt=0 tracks
  if (fabs(pxDF) < small) pxDF = pxDF<0?-small:small; // hack to avoid pt=0 tracks
  invpt = 1./invpt;  
    
  //omega
  double Bval = fieldmap.bFieldNominal();
  //  if (fabs(Bval) < small) {  // hack to avoid B=0 tracks (very far Z)
  //    Bval = small;
  //  }

//  pars(3) = -BField::cmTeslaToGeVc*Bval*sign/pt;
  pars(3) = invpt;
  pars(3) *= sign;
  pars(3) *= Bval;
  pars(3) *= -BField::cmTeslaToGeVc;

//  pars(5) = pzDF / pt;  //tandip
  pars(5) = pzDF;
  pars(5) *= invpt;

  DifNumber rho = 1./pars[2];
  phip=atan2(pyDF,pxDF);
  phip.cosAndSin(cosphip,sinphip);

//  pars(1) = (rho + yDF*cosphip - xDF*sinphip)/cos(gamma) - rho;//d0
  pars(1) = rho;
  pars(1) += yDF*cosphip;
  pars(1) -= xDF*sinphip;  // continued below ...

  gamma=atan((xDF*cosphip+yDF*sinphip)/ -pars(1)); // NOTE: do NOT use atan2 here
    
  pars(1) /= cos(gamma);
  pars(1) -= rho;

//  pars(2) = phip+gamma;  //phi0
  pars(2) = phip;
  pars(2) += gamma;
  //  pars(2).mod(0., Constants::twoPi);


//  pars(4) = zDF + rho*gamma*pars[4];  //z0
  pars(4) = pars[4];  // weird
  pars(4) *= gamma;
  pars(4) *= rho;
  pars(4) += zDF;


// Get error matrix for position and momentum

  static HepSymMatrix posandmomErr(6);
  static HepVector parsVec(5);

  int i;
  for (i = 1; i <= 3; ++i) {
    int j;
    for (j = 1; j <= i; ++j) {
      // with "fast" access, row should be >= col
      posandmomErr.fast(i,j) = pos.covMatrix().fast(i,j);
      posandmomErr.fast(i+3,j+3) = pmom.covMatrix().fast(i,j);
    }
    for (j = 1; j <= 3; ++j) {
      posandmomErr.fast(j+3,i) = cxp(i,j);
    }
  }
  for (i = 1; i <= 5; ++i) {
 // make the array of DifNums into a HepVector
 // (needed for TrkExchangePar init)
    parsVec(i) = pars(i).number();
  }
// Now calculate error on the helix pars--the real calculation

  return TrkExchangePar(parsVec,posandmomErr.similarity(pars.jacobian()) );
}
//----------------------------------------------------------------------
NeutParams TrkHelixUtils::lineFromMomErr(const BbrPointErr& pos,
              const BbrVectorErr& pmom,const HepMatrix& cxp, double sign, 
              const BField& fieldmap) {
//----------------------------------------------------------------------

  DifNumber xDF(pos.x(),1,6), yDF(pos.y(),2,6), zDF(pos.z(),3,6);
  DifNumber pxDF(pmom.x(), 4, 6);
  DifNumber pyDF(pmom.y(), 5, 6);
  DifNumber pzDF(pmom.z(), 6, 6);

  static DifArray pars(6,6);

  DifNumber pt = pxDF;
  pt *= pxDF;
  pt += pyDF*pyDF;

  if (pt < small) pt = small;               // hack to avoid pt=0 tracks
  if (fabs(pxDF) < small) pxDF = pxDF<0?-small:small; // hack to avoid pt=0 tracks

//  pars(3) = sqrt(pt*pt+pzDF*pzDF); // Magnitude of p 
  pars(3) = pzDF;
  pars(3) *= pzDF;
  pars(3) += pt;
  pars(3) = sqrt(pars(3));

  pt = sqrt(pt);
  DifNumber invpt = 1./pt;

  //Next lines modified by Eugenio Paoloni 19-XI-98

//  DifNumber pvxDF=pxDF/pt; // versor along pt x
  DifNumber pvxDF=pxDF;
  pvxDF *= invpt;

//  DifNumber pvyDF=pyDF/pt; // and y component
  DifNumber pvyDF=pyDF;
  pvyDF *= invpt;


//  pars(5) = pzDF / pt;  //tandip
  pars(5) = pzDF;
  pars(5) *= invpt;

  pars(2) = atan2(pyDF,pxDF);  //phi0

//  pars(1) = yDF*pvxDF - xDF*pvyDF;//d0
  pars(1) = yDF;
  pars(1) *= pvxDF;
  pars(1) -= xDF*pvyDF;

//  pars(4) = zDF -pars(5)*(xDF*pvxDF+yDF*pvyDF) ;  //z0
  pars(4) = xDF;
  pars(4) *= pvxDF;
  pars(4) += yDF*pvyDF;

// insert this in the middle for speed
//  pars(6) = (xDF*pvxDF+yDF*pvyDF)*pars(3)/pt;//s0
  pars(6) = pars(4);
  pars(6) *= pars(3);
  pars(6) *= invpt;


  pars(4) *= -pars(5);
  pars(4) += zDF;




// Get error matrix for position and momentum

  static HepSymMatrix posandmomErr(6);
  static HepVector parsVec(6);

  int i;
  for (i = 1; i <= 3; ++i) {
    int j;
    for (j = 1; j <= i; ++j) {
      // with "fast" access, row should be >= col
      posandmomErr.fast(i,j) = pos.covMatrix().fast(i,j);
      posandmomErr.fast(i+3,j+3) = pmom.covMatrix().fast(i,j);
    }
    for (j = 1; j <= 3; ++j) {
      posandmomErr.fast(j+3,i) = cxp(i,j);
    }
  }
  for (i = 1; i <= 6; ++i) {
 // make the array of DifNums into a HepVector
 // (needed for TrkExchangePar init)
    parsVec(i) = pars(i).number();
  }
// Now calculate error on the helix pars--the real calculation
  return NeutParams(parsVec,posandmomErr.similarity(pars.jacobian()) );
}

//------------------------------------------------------------------------
double TrkHelixUtils::fltToRad(const TrkExchangePar& hel, double rad) {
//------------------------------------------------------------------------
  double d0 = hel.d0();
  double omega = hel.omega();
  double tanDip = hel.tanDip();

  // GS To be able to use this with straight lines:  
//  if( fabs(omega) < 0.0001 ) omega=copysign(0.0001,omega);  // assume the pt= 1000 GeV  
  if( fabs(omega) < 0.0001 ) omega = (omega<0.0) ? -0.0001 : 0.0001 ;  
  
  double stuff =  ( rad*rad - d0*d0 ) / (1 + omega * d0);
  if (stuff <= 0.0) return 0.;
  if (omega == 0.) return sqrt(stuff);
  double sinAngle = 0.5 * omega * sqrt(stuff);
  double dist2d = 0;
  if (sinAngle < -1.0 || sinAngle > 1.0) {
    dist2d = Constants::pi / fabs(omega);
  } else {
    dist2d = 2. * asin(sinAngle) / omega;
  }
  return dist2d * sqrt(1. + tanDip*tanDip);
}
