//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: BbrPointErr.hh 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//      Endow a Point class with a covariance matrix
//      
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Forest Rouse
//      Victoria Novotny  
//      Ed Frank             University of Penn., efrank@upenn5.hep.upenn.edu
//
// History
//      14 Oct 96 Ed Frank   Simple mod to make unary - compile under
//                           native SunOS5 compiler.
//
// Copyright Information:
//	Copyright (C) 1996 U.C. Davis
//
//------------------------------------------------------------------------
#ifndef BBRPOINTERR_HH
#define BBRPOINTERR_HH

#include "BaBar/BaBar.hh"
#include "BbrGeom/BbrVectorErr.hh"
#include "BbrGeom/BbrError.hh"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Geometry/Translation.h"
#include "CLHEP/Geometry/Transformation.h"

class BbrPointErr : public HepPoint
{
public:
  enum PolarCoordinateIndex {   
    Rho = 0, 
    Theta = 1, 
    Phi = 2,
    NUM_PCOORDINATES = 3
  };
  enum CylindricalCoordinateIndex {   
    C_Rho = 0, 
    C_Zeta = 1, 
    C_Phi = 2,
    NUM_CCOORDINATES = 3
  };
  BbrPointErr() : HepPoint(), _covMatrix(3)		{}
  BbrPointErr(const HepPoint &p) : HepPoint(p),
  _covMatrix(3)						{}
  BbrPointErr(const HepPoint &p, const BbrError& covMat) : HepPoint(p),
  _covMatrix(3)						{ _covMatrix=covMat; }

  // new constructors for this class
  BbrPointErr(const BbrPointErr& v) : _covMatrix(3)	{*this = v;}

  BbrPointErr& operator=(const BbrPointErr& v)
    {
      if (this != &v) {
	HepPoint::operator=(v);
	_covMatrix = v.covMatrix();
      }
      return *this;
    }

  BbrPointErr operator - () {
      HepPoint t = *this;
      return BbrPointErr( -t, _covMatrix); // -covMatrix remains unaltered
  }

  BbrPointErr& operator -= (const BbrVectorErr& v){
      HepPoint::operator -= (v);
      _covMatrix += v.covMatrix();
      return *this;
  }

  BbrPointErr& operator += (const BbrVectorErr& v){
      HepPoint::operator += (v);
      _covMatrix += v.covMatrix();
      return *this;
  }

  BbrPointErr& transform(const HepTranslation& trans){
      HepPoint::transform(trans);
      return *this;
  }

  BbrPointErr& transform(const HepRotation& rot){
      HepPoint::transform(rot);
      _covMatrix = _covMatrix.similarity(rot);
      return *this;
  }

  BbrPointErr& transform(const HepTransformation& transf){
      HepPoint::transform(transf);
      _covMatrix = _covMatrix.similarity(transf.rot_mat());
      return *this;
  }

  // destructor MAY be needed later
  // virtual ~BbrPointErr() {};

  double determineChisq(const HepPoint& diffPoint) const
    {
      HepVector temp(3);
      temp[0] = diffPoint.x()-this->x();
      temp[1] = diffPoint.y()-this->y();
      temp[2] = diffPoint.z()-this->z();
      return _covMatrix.determineChisq(temp);
    }

  BbrError covRTPMatrix() const;
  // returns the covariance Matrix in spherical coordinate
  // use   PolarCoordinateIndex enum to get the components
  BbrError covRZPMatrix() const;
  // returns the covariance Matrix in cylindrical coordinate
  // use   CylindricalCoordinateIndex enum to get the components

  inline const BbrError& covMatrix() const		{ return _covMatrix; }
  inline void setCovMatrix(const BbrError& v)		{ _covMatrix = v; }

//  void printOn(ostream& out=cout) const;

private:
  
  BbrError _covMatrix;
};

BbrPointErr operator + (const BbrPointErr&, const BbrVectorErr&);

BbrPointErr operator + (const BbrVectorErr&, const BbrPointErr&);

BbrPointErr operator - (const BbrPointErr&, const BbrVectorErr&);

BbrVectorErr operator - (const BbrPointErr&, const BbrPointErr&);

std::ostream & operator<<(std::ostream & stream, const BbrPointErr & verr);
std::istream & operator>>(std::istream & stream, BbrPointErr & verr);


#endif





