// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetAlignElem.cc,v 1.21 2004/08/06 05:58:28 bartoldu Exp $
//
//  Description:  Define the alignment parameters for a detector element.  These
//  consist of a translation vector plus a set of alignment angles (sequential
//  rotations about x,y,z axes).  Access in the form of a HepTransform is
//  provided, but explicit application of the alignment must be coded according
//  to specific detector needs.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 12/5/96
//------------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/Constants.hh"
#include "BTrk/BaBar/ErrLog.hh"
#include "BTrk/DetectorModel/DetAlignElem.hh"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "BTrk/BbrGeom/AngleSets.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "BTrk/BbrGeom/Transformation.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include <string>
using std::endl;
using std::ostream;

//
//  Constructors
//
DetAlignElem::DetAlignElem(const Hep3Vector& dvec,const AlignAngles& rvec,
			   const char* name,int ielem):
  _ielem(ielem)
{
  _parvec[alndx] = dvec.x();
  _parvec[alndy] = dvec.y();
  _parvec[alndz] = dvec.z();
  _parvec[alnrx] = rvec.a;
  _parvec[alnry] = rvec.b;
  _parvec[alnrz] = rvec.c;
  std::string sname(name);
// we can only handle 19 characters in a DetAlignElem name
  if(sname.length() > NAMELENGTH){
    ErrMsg(error) << "DetAlignElem:Error, given name is > "
	 << NAMELENGTH << " characters, truncting" << endmsg;
   sname.resize(NAMELENGTH);
  }
  strcpy(_ename,sname.c_str());
  for(int icov=0;icov<21;icov++)
    _parcov[icov] = 0.0;
}

//
//  Construct from the parameter vector and covariance
//
DetAlignElem::DetAlignElem(const HepVector& pvec,const HepSymMatrix& pcov,
			   const char* name,int ielem):
  _ielem(ielem)
{
//
//  Check that the dimensions were all setup OK
//
  assert(6 == pvec.num_row());
//
//  Fill the data members
//
  for(int ipar=0;ipar<6;ipar++)
    _parvec[ipar] = pvec[ipar];
  setCovariance(pcov);
  strcpy(_ename,name);
}

//
//  construct from a HepTransformation (error set to 0)
//
DetAlignElem::DetAlignElem(const HepTransformation & trans,
			   const char* name,int ielem):
  _ielem(ielem)
{

  HepTranslation tl = trans.tra_vec();
  _parvec[0] = tl.x();
  _parvec[1] = tl.y();
  _parvec[2] = tl.z();

  HepRotation rt = trans.rot_mat();


  strcpy(_ename,name);
}

//
//  Copy constructor
//
DetAlignElem::DetAlignElem(const DetAlignElem& other):
  _ielem(other._ielem)
{
  for(int ipar=0;ipar<6;ipar++)
    _parvec[ipar] = other._parvec[ipar];
  for(int icov=0;icov<21;icov++)
    _parcov[icov] = other._parcov[icov];
 strcpy(_ename,other._ename);
}
//
//  Dummy constructors (fill data elements with 0s)
//
DetAlignElem::DetAlignElem(const char* name,int id):
  _ielem(id)
{
  for(int ipar=0;ipar<6;ipar++)
    _parvec[ipar] = 0.0;
  for(int icov=0;icov<21;icov++)
    _parcov[icov] = 0.0;
  strcpy(_ename,name);
}
DetAlignElem::DetAlignElem():
  _ielem(-1)
{
  for(int ipar=0;ipar<6;ipar++)
    _parvec[ipar] = 0.0;
  for(int icov=0;icov<21;icov++)
    _parcov[icov] = 0.0;
  strcpy(_ename,"Unknown");
}
//
//  Operators
//
DetAlignElem& DetAlignElem::operator = (const DetAlignElem& other){
  if(this != &other){
    _ielem = other._ielem;
    strcpy(_ename,other._ename);
    for(int ipar=0;ipar<6;ipar++)
      _parvec[ipar] = other._parvec[ipar];
    for(int icov=0;icov<21;icov++)
      _parcov[icov] = other._parcov[icov];
  }
  return *this;
}

bool
DetAlignElem::operator == (const DetAlignElem& other) const {
  return elementName() == other.elementName() &&
    elementNumber() == other.elementNumber();
}


//
//  Pretty printout
//
void DetAlignElem::print(ostream& out) const {
  out << "Alignment for element " << _ielem << _ename << endl;
  out << " Translation vector = " << displacement() << endl;
  out << " Rotation vector = " << "("<< _parvec[alnrx] << ", "
      << _parvec[alnry]   << ", " << _parvec[alnrz] <<")"<< endl;
}
//
//  Vector and convariance matrix of parameters, for algebraic manipulation
//
HepVector
DetAlignElem::parameterVector() const {
  HepVector pvec(6);
  for(int ipar=0;ipar<6;ipar++)
    pvec[ipar] = _parvec[ipar];
  return pvec;
}
//
HepSymMatrix
DetAlignElem::parameterCovariance() const {
  HepSymMatrix covar(6,0);
  int icov=0;
  for(int irow=0;irow<6;irow++){
    for(int icol=0;icol<=irow;icol++){
      covar.fast(irow+1,icol+1) = _parcov[icov++];
    }
  }
  assert(icov == 21);
  return covar;
}

Hep3Vector
DetAlignElem::displacement() const { return Hep3Vector(_parvec[alndx],
						       _parvec[alndy],_parvec[alndz]); }
AlignAngles
DetAlignElem::angles() const { return AlignAngles(_parvec[alnrx],
						  _parvec[alnry],_parvec[alnrz]); }
HepTransformation
DetAlignElem::transform() const { return HepTransformation(displacement(),
							   angles()); }

HepTransformation
DetAlignElem::inverseTransform() const { return transform().invert();}

std::string
DetAlignElem::elementName() const {return std::string(_ename); }

void
DetAlignElem::setCovariance(const HepSymMatrix& covar) {
  assert(6 == covar.num_row() &&
	 6 == covar.num_col() );
  int icov=0;
  for(int irow=0;irow<6;irow++){
    for(int icol=0;icol<=irow;icol++){
      _parcov[icov++] = covar.fast(irow+1,icol+1);
    }
  }
  assert(icov == 21);
}

// DetAlignElem multiplication
// attention: a *= b means a = a * b !!!
// the propagation of the coavriance matrix is not implemented !!!
DetAlignElem &
DetAlignElem::operator*=( const DetAlignElem & b ) {



  // translational part
  _parvec[alndx] += b._parvec[alndx];
  _parvec[alndy] += b._parvec[alndy];
  _parvec[alndz] += b._parvec[alndz];


  // rotational part
  HepRotation ar;
  HepRotation br;

  ar.rotateX( _parvec[alnrx] );
  ar.rotateY( _parvec[alnry] );
  ar.rotateZ( _parvec[alnrz] );

  br.rotateX( b._parvec[alnrx] );
  br.rotateY( b._parvec[alnry] );
  br.rotateZ( b._parvec[alnrz] );

  ar *= br;

  AlignAngles aa = alignAngles( ar );

  _parvec[alnrx] = aa.a;
  _parvec[alnry] = aa.b;
  _parvec[alnrz] = aa.c;

  return *this;
}

//
// private member functions
//

// this function returns the set of angles 'AlignAngles' from a given
// HepRotation - needed for the *= operator of DetAlignElem
// the problem of determining the AlignAngles from a HepRotation is unique if
// you require the sum of angle alpha and gamma ( around x and z respectively ) to be smaller than Pi
//
// AlignAngles( a, b, c ) defines a rotation around the axis' x, y, and z by angles a, b, and c in that order
AlignAngles
DetAlignElem::alignAngles( const HepRotation & r ) {

  float cosb2 = 1-r.zx()*r.zx();

  float sinapc = (r.xx()*r.zy()+r.yx()*r.zz())/cosb2;
  float cosapc = (r.xx()*r.zz()-r.yx()*r.zy())/cosb2;


  float apc = angle( sinapc, cosapc );

  float sinamc = (r.xx()*r.zy()-r.yx()*r.zz())/cosb2;
  float cosamc = (r.xx()*r.zz()+r.yx()*r.zy())/cosb2;

  float amc = angle( sinamc, cosamc );

  float a = (apc+amc)/2;
  float c = (apc-amc)/2;
  //  float b = angle( -r.zx(), r.zy()/sin(a) );
  float b = angle( -r.zx(), r.xx()/cos(c) );

  AlignAngles aa( a, b, c );

  return aa;
}

// angle from given sin and the cos of the angle
float
DetAlignElem::angle( float sina, float cosa ) {

  int signcos = 0;

  if( cosa>0 ) signcos = 1;
  if( cosa<0 ) signcos = -1;

  float a=asin(sina);
  if( signcos == 1 ) return a;
  if( sina > 0 ) return Constants::pi-a;
  return -(Constants::pi+a);
}

std::string
DetAlignElem::elName() const
{
  return elementName();
}
