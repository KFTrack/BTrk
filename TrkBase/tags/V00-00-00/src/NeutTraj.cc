#include "BaBar/BaBar.hh"
#include <assert.h>
#include <math.h>

#include "BaBar/Constants.hh"
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "TrkBase/NeutTraj.hh"
#include "TrkBase/TrkVisitor.hh"
#include "difAlgebra/DifNumber.hh"
#include "difAlgebra/DifPoint.hh"
#include "difAlgebra/DifVector.hh"
#include "BbrGeom/BbrAngle.hh"

NeutTraj::NeutTraj(const NeutParams& np, double lowlim,double hilim,
                   const HepPoint& refpoint) : 
  TrkSimpTraj(np.parameter(), np.covariance(), lowlim,hilim)
{
}

NeutTraj::NeutTraj( const NeutTraj& n )
  : TrkSimpTraj(n.parameters()->parameter(), n.parameters()->covariance(),
                n.lowRange(),n.hiRange())
{
}

//  Simple cloning function
NeutTraj*
NeutTraj::clone() const
{
  return new NeutTraj(*this);
}

NeutTraj&
NeutTraj::operator=(const NeutTraj& n)
{
  if( &n != this ){
    for(int iend=0;iend<2;iend++)
      flightrange[iend] = n.flightrange[iend];
    _dtparams = *n._np();
  }
  return *this;
}

//-----------------------------------------------------------------
NeutTraj::~NeutTraj() {
}

//-----------------------------------------------------------------
double 
NeutTraj::x( const double& f ) const {
//-----------------------------------------------------------------
  return cosDip()*f*cos(phi0()) - d0()*sin(phi0());
}

//-----------------------------------------------------------------
double 
NeutTraj::y( const double& f ) const {
//-----------------------------------------------------------------
  return cosDip()*f*sin(phi0()) + d0()*cos(phi0());
}

//-----------------------------------------------------------------
double 
NeutTraj::z(const double& f) const {
//-----------------------------------------------------------------
  return (z0() + f*sinDip());
}

//-----------------------------------------------------------------
HepPoint 
NeutTraj::position( double f) const {
//-----------------------------------------------------------------
  double cdip = cosDip();
  double sphi0 = sin(phi0());
  double cphi0 = cos(phi0());

  return HepPoint(cdip * f * cphi0 - d0()*sphi0, 
		  cdip * f * sphi0 + d0()*cphi0, 
		  z0() + f * tanDip() * cdip);
}

//-----------------------------------------------------------------
Hep3Vector   
NeutTraj::direction( double f) const {
//-----------------------------------------------------------------
  return Hep3Vector ( cos(phi0())*cosDip(), sin(phi0())*cosDip(),
		               sinDip() );
}
//-----------------------------------------------------------------
Hep3Vector 
NeutTraj::delDirect( double fltLen ) const{
//-----------------------------------------------------------------
  return Hep3Vector(0.0, 0.0, 0.0);
}

//-----------------------------------------------------------------
double 
NeutTraj::distTo1stError(double , double tol, int) const {
//-----------------------------------------------------------------
  return 9999.;
}

//-----------------------------------------------------------------
double
NeutTraj::distTo2ndError(double , double tol, int) const {
//-----------------------------------------------------------------
  return 9999.;
}

void
NeutTraj::getInfo(double fltLen, HepPoint& pos, Hep3Vector& dir,
                  Hep3Vector& delDir) const
{
  // This could be made much more efficient!!!!!!
  pos = position(fltLen);
  dir = direction(fltLen);
  delDir = delDirect(fltLen);
}

void
NeutTraj::getInfo( double fltLen, HepPoint& pos, Hep3Vector& dir ) const
{
  //  This could be made much more efficient!!!!!
  pos = position(fltLen);
  dir = direction(fltLen);
  return;
}

void
NeutTraj::getDFInfo(double flt, DifPoint& pos, DifVector& dir,
                    DifVector& delDir) const
{
  //Provides difNum version of information for calculation of derivatives.
  //  Some of this duplicate code could be reduced if we had member 
  //  function templates.  Maybe.

  // Create difNumber versions of parameters
  //enum index (phi0(), etc) is from HelixParams.hh
  DifNumber phi0Df(phi0(), NeutParams::_phi0+1, NeutParams::_nneutprm);
  DifNumber d0Df(d0(), NeutParams::_d0+1, NeutParams::_nneutprm);
  DifNumber z0Df(z0(), NeutParams::_z0+1, NeutParams::_nneutprm);
  DifNumber tanDipDf(tanDip(), NeutParams::_tanDip+1, NeutParams::_nneutprm);
  DifNumber pDf(p(), NeutParams::_p+1, NeutParams::_nneutprm);
  // RF 03/25/99
  // the flight lenght has an error, which is the error along the trajectory 
  // in the creation point
  DifNumber s0Df(flt, NeutParams::_s0+1, NeutParams::_nneutprm);
  phi0Df.setIndepPar( parameters() );
  d0Df.setIndepPar( parameters() );
  z0Df.setIndepPar( parameters() );
  tanDipDf.setIndepPar( parameters() );
  pDf.setIndepPar( parameters() );
  s0Df.setIndepPar( parameters() );
  DifNumber dipDf = atan(tanDipDf);

  DifNumber sDip, cDip;
  dipDf.cosAndSin(cDip, sDip);
  DifNumber sinPhi0, cosPhi0;
  phi0Df.cosAndSin(cosPhi0, sinPhi0);

  DifNumber x =  cDip*s0Df*cosPhi0 - d0Df*sinPhi0;
  DifNumber y =  cDip*s0Df*sinPhi0 + d0Df*cosPhi0;
  DifNumber z =  z0Df + s0Df*sDip ;
  pos =  DifPoint(x, y, z);
  dir = DifVector( cosPhi0*cDip, sinPhi0*cDip, sDip );
  delDir = DifVector(0., 0., 0.);
}

HepMatrix
NeutTraj::derivDeflect(double fltlen,deflectDirection idirect) const
{
//
//  This function computes the column matrix of derivatives for the change
//  in parameters for a change in the direction of a track at a point along
//  its flight, holding the momentum and position constant.  The effects for
//  changes in 2 perpendicular directions (theta1 = dip and 
//  theta2 = phi*cos(dip)) can sometimes be added, as scattering in these
//  are uncorrelated.
//
  HepMatrix ddflct(NeutParams::_nneutprm,1); 
//
//  Compute some common things
//
  double tand = tanDip();
  double cosd = cosDip();
//
//  Go through the parameters
//
  switch (idirect) {
  case theta1:
    ddflct[NeutParams::_p][0] = 0.;
    ddflct[NeutParams::_tanDip][0] = 1.0/sqr(cosd);
    ddflct[NeutParams::_d0][0] = 0.;
    ddflct[NeutParams::_phi0][0] =  0.;
    ddflct[NeutParams::_s0][0] =  0.;
    ddflct[NeutParams::_z0][0] = translen(fltlen) * (-1. - sqr(tand));
    break;
  case theta2:
    ddflct[NeutParams::_p][0] = 0;
    ddflct[NeutParams::_tanDip][0] = 0;
    ddflct[NeutParams::_d0][0] = -translen(fltlen)/cosd;
    ddflct[NeutParams::_phi0][0] = 1./cosd;
    ddflct[NeutParams::_z0][0] = 0.;
    ddflct[NeutParams::_s0][0] =  0.;
    break;
  }
  return ddflct;
}

HepMatrix
NeutTraj::derivDisplace(double fltlen,deflectDirection idirect) const
{
//
//  This function computes the column matrix of derivatives for the change
//  in parameters for a change in the transvers position of a track at a point along
//  its flight, holding the momentum and direction constant.  The effects for
//  changes in 2 perpendicular directions (theta1 = dip and 
//  theta2 = phi*cos(dip)) are uncorrelated.  The sign convention has been
//  chosen so that correlated scattering and displacement effects may be added
//
  HepMatrix ddflct(NeutParams::_nneutprm,1); 
//
//  Compute some common things
//
  double cosd = cosDip();
//
//  Go through the parameters
//
  switch (idirect) {
  case theta1:
    ddflct[NeutParams::_p][0] = 0.;
    ddflct[NeutParams::_tanDip][0] = 0.0;
    ddflct[NeutParams::_d0][0] = 0.;
    ddflct[NeutParams::_phi0][0] =  0.;
    ddflct[NeutParams::_s0][0] =  0.;
    ddflct[NeutParams::_z0][0] = 1.0/cosd;
    break;
  case theta2:
    ddflct[NeutParams::_p][0] = 0;
    ddflct[NeutParams::_tanDip][0] = 0;
    ddflct[NeutParams::_d0][0] = 1.0;
    ddflct[NeutParams::_phi0][0] = 0.0;
    ddflct[NeutParams::_z0][0] = 0.;
    ddflct[NeutParams::_s0][0] =  0.;
    break;
  }
  return ddflct;
}

HepMatrix
NeutTraj::derivPFract(double fltlen) const
{
//  This function computes the column matrix of derivatives for the change
//  in parameters from a (fractional) change in the track momentum,
//  holding the direction and position constant.  The momentum change can
//  come from energy loss or bfield inhomogeneities.
//
//  For a helix, dp/P = -domega/omega,
//  dParam/d(domega/omega) = -omega*dParam/ddomega
//
  HepMatrix dmomfrac(NeutParams::_nneutprm,1);

//  Go through the parameters
  dmomfrac[NeutParams::_p][0] = p(); // momentum
  dmomfrac[NeutParams::_tanDip][0] = 0.0; // tanDip
  dmomfrac[NeutParams::_d0][0] = 0.0; // d0
  dmomfrac[NeutParams::_phi0][0] = 0.0; // phi0
  dmomfrac[NeutParams::_z0][0] = 0.0; // z0
  dmomfrac[NeutParams::_s0][0] = 0.0; // s0
  return dmomfrac;
}


double
NeutTraj::phi0() const
{
  return BbrAngle(_np()->phi0()).rad();
}

void
NeutTraj::paramFunc(const HepPoint& oldpoint,const HepPoint& newpoint,
                    const HepVector& oldvec,const HepSymMatrix& oldcov,
                    HepVector& newvec,HepSymMatrix& newcov,
                    double fltlen)
{

// copy the input parameter vector, in case the input and output are the same
  HepVector parvec(oldvec);
// start with the input: momentum and tandip don't change
  newvec = parvec;
//
  double delx = newpoint.x()-oldpoint.x();
  double dely = newpoint.y()-oldpoint.y();
  double delz = newpoint.z()-oldpoint.z();
//  
  double cos0 = cos(parvec[NeutParams::_phi0]);
  double sin0 = sin(parvec[NeutParams::_phi0]);
  double perp = delx*sin0-dely*cos0;
  double para = delx*cos0+dely*sin0;
  double tand = parvec[NeutParams::_tanDip];
// delta
// assume delta, newdelta have the same sign
// d0
  newvec[NeutParams::_d0] = parvec[NeutParams::_d0] + perp;
// phi0; check that we don't get the wrong wrapping
  newvec[NeutParams::_phi0] = parvec[NeutParams::_phi0];
//z0   
  newvec[NeutParams::_z0] += 
    parvec[NeutParams::_tanDip]*(delx/cos0)*(1.+(sin0/cos0)) - delz;
// now covariance: first, compute the rotation matrix
  HepMatrix covrot(NeutParams::_nneutprm,NeutParams::_nneutprm,0); // start with 0: lots of terms are zero
//
// momentum is diagonal
  covrot[NeutParams::_p][NeutParams::_p] = 1.0;
// tandip is diagonal
  covrot[NeutParams::_tanDip][NeutParams::_tanDip] = 1.0;
// d0
  covrot[ NeutParams::_d0][ NeutParams::_d0] = 1.;
  covrot[ NeutParams::_d0][ NeutParams::_phi0] = para;
// phi0
  covrot[ NeutParams::_phi0][ NeutParams::_p] = para;
  covrot[ NeutParams::_phi0][ NeutParams::_phi0] = 1.;
// z0
  covrot[ NeutParams::_z0][ NeutParams::_phi0] = tand*-2.*perp;
  covrot[ NeutParams::_z0][ NeutParams::_tanDip] = 
                                      (delx/cos0)*(1.+(sin0/cos0));
  covrot[ NeutParams::_z0][ NeutParams::_z0] = 1.0;  
  covrot[ NeutParams::_s0][ NeutParams::_s0] = 1.0;  
//  
//  Apply the rotation
  newcov = oldcov.similarity(covrot);
// done
}

void
NeutTraj::invertParams(TrkParams* newparams, std::vector<bool>& flags) const
{
  assert(1==0);

//  // Inverts parameters and returns true if the parameter inversion
//  // requires a change in sign of elements in the covariance matrix
//
//  for (unsigned iparam = 0; iparam < NeutParams::_nneutprm; iparam++) {
//    switch ( iparam ) {
//    case NeutParams::_d0:  // changes sign
//    case NeutParams::_tanDip:  // changes sign
//      newparams->parameter()[iparam] *= -1.0;
//      flags[iparam] = true;
//      break;
//    case NeutParams::_phi0:  // changes by pi, but covariance
//                             // matrix shouldn't changed
//      newparams->parameter()[iparam] =
//	BbrAngle(newparams->parameter()[iparam] + Constants::pi);
//      flags[iparam] = false;
//      break;
//    case NeutParams::_p:  // no change sign
//    case NeutParams::_z0:  // no change
//      flags[iparam] = false;
//    }
//  }
//  return;
}

void
NeutTraj::visitAccept(TrkVisitor* vis) const
{
// Visitor access--just use the Visitor class member function
  vis->trkVisitNeutTraj(this);
}
