// $Id: HelixTraj.cc,v 1.7 2004/09/30 13:39:06 hulsberg Exp $

#include "BaBar/BaBar.hh"
#include <assert.h>
#include <math.h>
#include <limits.h>

#include "BaBar/Constants.hh"
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkVisitor.hh"
#include "difAlgebra/DifNumber.hh"
#include "difAlgebra/DifPoint.hh"
#include "difAlgebra/DifVector.hh"
#include "BbrGeom/BbrAngle.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "ErrLogger/ErrLog.hh"
using std::endl;
using std::ostream;
//  Fix for some machines
//
#ifndef M_2PI
#define M_2PI 2*M_PI
#endif

HelixTraj::HelixTraj(const HepVector& pvec, const HepSymMatrix& pcov,
                     double lowlim, double hilim, const HepPoint& refpoint) :
  TrkSimpTraj(pvec, pcov, lowlim,hilim,refpoint)
{
  //  Make sure the dimensions of the input matrix and vector are correct

  if( pvec.num_row() != NHLXPRM || pcov.num_row() != NHLXPRM ){
    ErrMsg(fatal) 
      << "HelixTraj: incorrect constructor vector/matrix dimension" << endmsg;
  }

  if (omega() == 0.0) parameters()->parameter()[omegaIndex] = 1.e-9;
}


HelixTraj::HelixTraj(const TrkExchangePar& inpar,
                     double lowlim, double hilim, const HepPoint& refpoint) :
  TrkSimpTraj(inpar.params(), inpar.covariance(), lowlim,hilim,refpoint)
{
  if (omega() == 0.0) parameters()->parameter()[omegaIndex] = 1.e-9;
}

HelixTraj::HelixTraj(const TrkParams& inpar,
                     double lowlim, double hilim, const HepPoint& refpoint) :
  TrkSimpTraj(inpar, lowlim,hilim,refpoint)
{
  assert(inpar.parameter().num_row()==NHLXPRM);
  if (omega() == 0.0) parameters()->parameter()[omegaIndex] = 1.e-9;
}

HelixTraj::HelixTraj( const HelixTraj& h )
  : TrkSimpTraj(h.parameters()->parameter(), h.parameters()->covariance(), 
                h.lowRange(),h.hiRange(),h.referencePoint())
{
}

HelixTraj*
HelixTraj::clone() const
{
  return new HelixTraj(*this);
}

HelixTraj&
HelixTraj::operator=(const HelixTraj& h)
{
  if( &h != this ){
    Trajectory::operator=(h);
    _dtparams = *(h.parameters());
    _refpoint = h._refpoint;
  }
  return *this;
}

HelixTraj::~HelixTraj()
{
}

double
HelixTraj::z(const double& f) const
{
  return z0() + f*sinDip() + referencePoint().z();
}

HepPoint
HelixTraj::position( double f) const
{
  double cDip = cosDip();
  double sDip = tanDip() * cDip;
  double phi00 = parameters()->parameter()[phi0Index];  // Don't normalize
  double ang = phi00 + cDip*f*omega();
  double cang = cos(ang);
  double sang = sin(ang);
  double sphi0 = sin(phi00);
  double cphi0 = cos(phi00);

  return HepPoint((sang - sphi0)/omega() - d0()*sphi0+referencePoint().x(),
                 -(cang - cphi0)/omega() + d0()*cphi0+referencePoint().y(),
                 z0() + f*sDip                       +referencePoint().z());
}

Hep3Vector
HelixTraj::direction( double f) const
{
  // Angle formed by tangent vector after
  // being rotated 'arclength' around orbit.
  double alpha = angle( f );
  // Construct 3-D tangent vector of unit magnitude.
  double cDip = cosDip();
  return Hep3Vector ( cos(alpha)*cDip,
                      sin(alpha)*cDip,
                      cDip*tanDip() );
}

Hep3Vector
HelixTraj::delDirect( double fltLen ) const
{
  double ang = angle(fltLen);
  double cDip = cosDip();
  double delX = -omega() * cDip * cDip * sin(ang);
  double delY =  omega() * cDip * cDip * cos(ang);
  return Hep3Vector(delX, delY, 0.0);
}

double
HelixTraj::distTo1stError(double s, double tol, int pathDir) const
{
  return sqrt(2.*tol/fabs(omega())*(1.+pow(tanDip(),2)));
}

double
HelixTraj::distTo2ndError(double s, double tol, int pathDir) const
{
  return sqrt(1.+pow(tanDip(),2))*cbrt(6.*tol/pow(omega(),2));
}

void
HelixTraj::getInfo(double fltLen, HepPoint& pos, Hep3Vector& dir,
                   Hep3Vector& delDir) const
{
  //  double ang = angle(fltLen);
  double cDip = cosDip();
  double sDip = tanDip() * cDip;
  double phi00 = parameters()->parameter()[phi0Index];  // Don't normalize
  double ang = phi00 + cDip*fltLen*omega(); 
  double cang = cos(ang);
  double sang = sin(ang);
  double sphi0 = sin(phi00);
  double cphi0 = cos(phi00);

  double xt = (sang - sphi0)/omega() - d0()*sphi0 +
    referencePoint().x();
  double yt = -(cang - cphi0)/omega() + d0()*cphi0 +
    referencePoint().y();
  double zt = z0() + fltLen*sDip + referencePoint().z();
  pos.setX(xt);
  pos.setY(yt);
  pos.setZ(zt);

  dir.setX(cang * cDip);
  dir.setY(sang * cDip);
  dir.setZ(sDip);

  double delX = -omega() * cDip * cDip * sang;
  double delY =  omega() * cDip * cDip * cang;
  delDir.setX(delX);
  delDir.setY(delY);
  delDir.setZ(0.0);
}

void
HelixTraj::getInfo( double fltLen, HepPoint& pos, Hep3Vector& dir ) const
{
  //  double ang = angle(fltLen);
  double cDip = cosDip();
  double sDip = tanDip() * cDip;
  double phi00 = parameters()->parameter()[phi0Index];  // Don't normalize
  double ang = phi00 + cDip*fltLen*omega(); 
  double cang = cos(ang);
  double sang = sin(ang);
  double sphi0 = sin(phi00);
  double cphi0 = cos(phi00);

  double xt = (sang - sphi0)/omega() - d0()*sphi0 +
    referencePoint().x();
  double yt = -(cang - cphi0)/omega() + d0()*cphi0 +
    referencePoint().y();
  double zt = z0() + fltLen*sDip + referencePoint().z();
  pos.setX(xt);
  pos.setY(yt);
  pos.setZ(zt);

  dir.setX(cang * cDip);
  dir.setY(sang * cDip);
  dir.setZ(sDip);
}

void HelixTraj::getDFInfo2(double flt, DifPoint& pos, DifVector& dir) const
{
  //Provides difNum version of information for calculation of derivatives.
  //  All arithmetic operations have been replaced by +=, etc. versions 
  //  for speed.

  // Create difNumber versions of parameters
  DifNumber phi0Df(phi0(), phi0Index+1, NHLXPRM);
  phi0Df.setIndepPar( parameters() );
  DifNumber d0Df(d0(), d0Index+1, NHLXPRM);
  d0Df.setIndepPar( parameters() );
  DifNumber z0Df(z0(), z0Index+1, NHLXPRM);
  z0Df.setIndepPar( parameters() );
  DifNumber tanDipDf(tanDip(), tanDipIndex+1, NHLXPRM);
  tanDipDf.setIndepPar( parameters() );
  DifNumber omegaDf(omega(), omegaIndex+1, NHLXPRM);
  omegaDf.setIndepPar( parameters() );

  DifNumber dipDf = atan(tanDipDf);

  static DifNumber cDip;
  dipDf.cosAndSin(cDip, dir.z);
  static DifNumber sinPhi0, cosPhi0;
  phi0Df.cosAndSin(cosPhi0, sinPhi0);

  bool lref = (referencePoint().x() != 0. || referencePoint().y() != 0. ||
               referencePoint().z() != 0.);

  DifNumber alphaDf = cDip;
  alphaDf *= omegaDf;
  alphaDf *= flt;
  alphaDf += phi0Df;

  // This is not the prettiest line imaginable for this operation:
  alphaDf.mod(-Constants::pi, Constants::pi);
  //  DifNumber sinAlpha, cosAlpha;
  alphaDf.cosAndSin(dir.x, dir.y);

  //  DifNumber x =   (sinAlpha - sinPhi0) / omegaDf - d0Df * sinPhi0 + px;
  //  DifNumber y =  -(cosAlpha - cosPhi0) / omegaDf + d0Df * cosPhi0 + py;

  pos.x = dir.y;
  pos.x -= sinPhi0;
  pos.x /= omegaDf;
  DifNumber temp = d0Df;
  temp *= sinPhi0;
  pos.x -= temp;

  pos.y = cosPhi0;
  pos.y -= dir.x;
  pos.y /= omegaDf;
  temp = d0Df;
  temp *= cosPhi0;
  pos.y += temp;

  pos.z = flt;
  pos.z *= dir.z;
  pos.z += z0Df;

  if (lref) {
    DifNumber px(referencePoint().x());
    DifNumber py(referencePoint().y());
    DifNumber pz(referencePoint().z());
    pos.x += px;
    pos.y += py;
    pos.z += pz;
  }

  dir.x *= cDip;
  dir.y *= cDip;
}

void
HelixTraj::getDFInfo(double flt, DifPoint& pos, DifVector& dir,
                     DifVector& delDir) const
{
  //Provides difNum version of information for calculation of derivatives.
  //  All arithmetic operations have been replaced by +=, etc. versions 
  //  for speed.

  // Create difNumber versions of parameters
  DifNumber phi0Df(phi0(), phi0Index+1, NHLXPRM);
  DifNumber d0Df(d0(), d0Index+1, NHLXPRM);
  DifNumber z0Df(z0(), z0Index+1, NHLXPRM);
  DifNumber tanDipDf(tanDip(), tanDipIndex+1, NHLXPRM);
  DifNumber omegaDf(omega(), omegaIndex+1, NHLXPRM);
  phi0Df.setIndepPar( parameters() );
  d0Df.setIndepPar( parameters() );
  z0Df.setIndepPar( parameters() );
  tanDipDf.setIndepPar( parameters() );
  omegaDf.setIndepPar( parameters() );
  DifNumber dipDf = atan(tanDipDf);

  static DifNumber cDip;
  dipDf.cosAndSin(cDip, dir.z);
  static DifNumber sinPhi0, cosPhi0;
  phi0Df.cosAndSin(cosPhi0, sinPhi0);

  bool lref = (referencePoint().x() != 0. || referencePoint().y() != 0. ||
	       referencePoint().z() != 0.);

  DifNumber alphaDf = cDip;
  alphaDf *= omegaDf;
  alphaDf *= flt;
  alphaDf += phi0Df;

  // This is not the prettiest line imaginable for this operation:
  alphaDf.mod(-Constants::pi, Constants::pi);
  //  DifNumber sinAlpha, cosAlpha;
  alphaDf.cosAndSin(dir.x, dir.y);

  //  DifNumber x =   (sinAlpha - sinPhi0) / omegaDf - d0Df * sinPhi0 + px;
  //  DifNumber y =  -(cosAlpha - cosPhi0) / omegaDf + d0Df * cosPhi0 + py;

  pos.x = dir.y;
  pos.x -= sinPhi0;
  pos.x /= omegaDf;
  DifNumber temp = d0Df;
  temp *= sinPhi0;
  pos.x -= temp;

  pos.y = cosPhi0;
  pos.y -= dir.x;
  pos.y /= omegaDf;
  temp = d0Df;
  temp *= cosPhi0;
  pos.y += temp;

  pos.z = flt;
  pos.z *= dir.z;
  pos.z += z0Df;

  if (lref) {
    DifNumber px(referencePoint().x());
    DifNumber py(referencePoint().y());
    DifNumber pz(referencePoint().z());
    pos.x += px;
    pos.y += py;
    pos.z += pz;
  }

  delDir.x = -omegaDf;
  delDir.x *= cDip;
  delDir.x *= cDip;
  delDir.x *= dir.y;

  delDir.y =  omegaDf;
  delDir.y *= cDip;
  delDir.y *= cDip;
  delDir.y *= dir.x;

  delDir.z = 0.;

  dir.x *= cDip;
  dir.y *= cDip;
}

HepMatrix
HelixTraj::derivDeflect(double fltlen,deflectDirection idirect) const
{
//
//  This function computes the column matrix of derrivatives for the change
//  in parameters for a change in the direction of a track at a point along
//  its flight, holding the momentum and position constant.  The effects for
//  changes in 2 perpendicular directions (theta1 = dip and
//  theta2 = phi*cos(dip)) can sometimes be added, as scattering in these
//  are uncorrelated.
//
  HepMatrix ddflct(NHLXPRM,1);
//
//  Compute some common things
//
  double omeg = omega();
  double tand = tanDip();
  double arcl = arc(fltlen);
  double dx = cos(arcl);
  double dy = sin(arcl);
  double cosd = cosDip();
  double darc = omeg*d0();
//
//  Go through the parameters
//
  switch (idirect) {
  case theta1:
    ddflct(omegaIndex+1,1) = omeg*tand;
    ddflct(tanDipIndex+1,1) = 1.0/pow(cosd,2);
    ddflct(d0Index+1,1) = (1-dx)*tand/omeg;
    ddflct(phi0Index+1,1) =  -dy*tand/(1+darc);
    ddflct(z0Index+1,1) = - translen(fltlen) - pow(tand,2)*dy/(omeg*(1+darc));
    break;
  case theta2:
    ddflct(omegaIndex+1,1) = 0;
    ddflct(tanDipIndex+1,1) = 0;
    ddflct(d0Index+1,1) = -dy/(cosd*omeg);
    ddflct(phi0Index+1,1) = dx/(cosd*(1+darc));
    ddflct(z0Index+1,1) = -tand*(1- dx/(1+darc))/(cosd*omeg);
    break;
  }

  return ddflct;
}


HepMatrix
HelixTraj::derivDisplace(double fltlen,deflectDirection idirect) const
{
//
//  This function computes the column matrix of derrivatives for the change
//  in parameters for a change in the position of a track at a point along
//  its flight, holding the momentum and direction constant.  The effects for
//  changes in 2 perpendicular directions 'theta1' = (-sin(l)cos(p),-sin(l)sin(p),cos(l)) and
//  'theta2' = (-sin(p),cos(p),0).  These are by definition orthogonal and uncorrelated.
//  these displacements are correlated with the angular change above
//
  HepMatrix ddflct(NHLXPRM,1);
//
//  Compute some common things
//
  double omeg = omega();
  double tand = tanDip();
  double arcl = arc(fltlen);
  double dx = cos(arcl);
  double dy = sin(arcl);
  double cosd = cosDip();
  double sind = sinDip();
  double darc_1 = 1.0+omeg*d0();
//
//  Go through the parameters
//
  switch (idirect) {
  case theta1:
    ddflct(omegaIndex+1,1) = 0.0;
    ddflct(tanDipIndex+1,1) = 0.0;
    ddflct(d0Index+1,1) = -sind*dy;
    ddflct(phi0Index+1,1) = sind*dx*omeg/darc_1;
    ddflct(z0Index+1,1) = sind*tand*dx/darc_1 + cosd;
    break;
  case theta2:
    ddflct(omegaIndex+1,1) = 0;
    ddflct(tanDipIndex+1,1) = 0;
    ddflct(d0Index+1,1) = dx;
    ddflct(phi0Index+1,1) = dy*omeg/darc_1;
    ddflct(z0Index+1,1) = tand*dy/darc_1;
    break;
  }

  return ddflct;
}


HepMatrix
HelixTraj::derivPFract(double fltlen) const
{
//
//  This function computes the column matrix of derrivatives for the change
//  in parameters from a (fractional) change in the track momentum,
//  holding the direction and position constant.  The momentum change can
//  come from energy loss or bfield inhomogeneities.
//
//  For a helix, dp/P = -domega/omega,
//  dParam/d(domega/omega) = -omega*dParam/ddomega
//
  HepMatrix dmomfrac(NHLXPRM,1);
//
//  Compute some common things

  double omeg = omega();
  double tand = tanDip();
  double tranl = translen(fltlen);
  double arcl = tranl*omeg;
  double dx = cos(arcl);
  double dy = sin(arcl);
  double darc = omeg*d0();

//  Go through the parameters
// omega
  dmomfrac(omegaIndex+1,1) = -omeg;
// tanDip
  dmomfrac(tanDipIndex+1,1) = 0.0;
// d0
  dmomfrac(d0Index+1,1) = -(1-dx)/omeg;
// phi0
  dmomfrac(phi0Index+1,1) = dy/(1+darc);
// z0
  dmomfrac(z0Index+1,1) = -tand*(tranl-dy/((1+darc)*omeg));
//
  return dmomfrac;
}

double
HelixTraj::curvature(double ) const
{
//  Compute the curvature as the magnitude of the 2nd derivative
//  of the position function with respect to the 3-d flight distance
//
  double cosd = cosDip();
  return pow(cosd,2)*fabs(omega());
}

double
HelixTraj::phi0() const
{
  return BbrAngle(parameters()->parameter()[phi0Index]).rad();
}

void
HelixTraj::paramFunc(const HepPoint& oldpoint,const HepPoint& newpoint,
                     const HepVector& oldvec,const HepSymMatrix& oldcov,
                     HepVector& newvec,HepSymMatrix& newcov,
                     double fltlen)
{
// copy the input parameter vector, in case the input and output are the same
  HepVector parvec(oldvec);
// start with the input: omega and tandip don't change
  newvec = parvec;
//
  double delx = newpoint.x()-oldpoint.x();
  double dely = newpoint.y()-oldpoint.y();
  double delz = newpoint.z()-oldpoint.z();
//
  double rad = 1./parvec[omegaIndex];
  double rad2 = rad*rad;
  double delta = rad + parvec[d0Index];
  double cos0 = cos(parvec[phi0Index]);
  double sin0 = sin(parvec[phi0Index]);
  double perp = delx*sin0-dely*cos0;
  double para = delx*cos0+dely*sin0;
  double tand = parvec[tanDipIndex];
  double oldphi  = parvec[phi0Index] +
    fltlen*parvec[omegaIndex]/sqrt(1.+tand*tand);
// delta
  double newdelta2 = delta*delta + delx*delx + dely*dely +
    2.0*delta*perp;
// assume delta, newdelta have the same sign
  double newdelta = delta>0 ? sqrt(newdelta2) : -sqrt(newdelta2);
  double invdelta = 1.0/newdelta;
  double invdelta2 = 1.0/newdelta2;
// d0
  newvec[d0Index] = newdelta - rad;
// phi0; check that we don't get the wrong wrapping. Atan2 has 2Pi ambiguity, not pi
  double newphi = atan2(sin0+delx/delta,cos0-dely/delta);
  while(fabs(newphi - oldphi)>M_PI)
    if(newphi > oldphi)
      newphi -= M_2PI;
    else
      newphi += M_2PI;
  newvec[phi0Index] = newphi;
  double delphi = newphi-parvec[phi0Index];
//z0
  newvec[z0Index] += tand*rad*(delphi) - delz;
// now covariance: first, compute the rotation matrix
// start with 0: lots of terms are zero
  static HepMatrix covrot(NHLXPRM,NHLXPRM,0);
//
// omega is diagonal
  covrot(omegaIndex+1,omegaIndex+1) = 1.0;
// tandip is diagonal
  covrot(tanDipIndex+1,tanDipIndex+1) = 1.0;
// d0
  covrot(d0Index+1,omegaIndex+1) = rad2*(1.0 - invdelta*(delta + perp));
  covrot(d0Index+1,d0Index+1) = invdelta*(delta + perp);
  covrot(d0Index+1,phi0Index+1) = delta*para*invdelta;
// phi0
  covrot(phi0Index+1,omegaIndex+1) = rad2*para*invdelta2;
  covrot(phi0Index+1,d0Index+1) = -para*invdelta2;
  covrot(phi0Index+1,phi0Index+1) = delta*(delta + perp)*invdelta2;
// z0
  covrot(z0Index+1,omegaIndex+1) = tand*
    (rad*covrot(phi0Index+1,omegaIndex+1) - rad2*delphi);
  covrot(z0Index+1,d0Index+1) = tand*rad*covrot(phi0Index+1,d0Index+1);
  covrot(z0Index+1,phi0Index+1) = 
    tand*rad*(covrot(phi0Index+1,phi0Index+1) - 1.0);
  covrot(z0Index+1,tanDipIndex+1) = rad*delphi;
  covrot(z0Index+1,z0Index+1) = 1.0;
//
//  Apply the rotation
  newcov = oldcov.similarity(covrot);
// done
}

void
HelixTraj::visitAccept(TrkVisitor* vis) const
{
// Visitor access--just use the Visitor class member function
  vis->trkVisitHelixTraj(this);
}

void
HelixTraj::invertParams(TrkParams* params, std::vector<bool>& flags) const
{
  // Inverts parameters and returns true if the parameter inversion
  // requires a change in sign of elements in the covariance matrix

  for (unsigned iparam = 0; iparam < NHLXPRM; iparam++) {
    switch ( iparam ) {
    case d0Index:  // changes sign
    case omegaIndex:  // changes sign
    case tanDipIndex:  // changes sign
      params->parameter()[iparam] *= -1.0;
      flags[iparam] = true;
      break;
    case phi0Index:  // changes by pi, but covariance matrix shouldn't change
      params->parameter()[iparam] =
        BbrAngle(params->parameter()[iparam] + Constants::pi);
      flags[iparam] = false;
      break;
    case z0Index:  // nochange
      flags[iparam] = false;
    }
  }
  return;
}

double
HelixTraj::angle(const double& f) const
{
  return BbrAngle(phi0() + arc(f));
}

void
HelixTraj::printAll(ostream& os) const
{
  os  << "HelixTraj with range "
      << lowRange() <<" to " << hiRange() << " and parameters " << endl
      << "d0= " << d0() << " phi0= "
      << phi0() << " omega= "
      << omega() << " z0 = "
      << z0() << " tanDip= "
      << tanDip() << endl;
}

void
HelixTraj::print(ostream& os) const
{
  Trajectory::print(os << "HelixTraj" );
}
