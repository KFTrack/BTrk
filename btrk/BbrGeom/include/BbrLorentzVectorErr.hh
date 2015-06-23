//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: BbrLorentzVectorErr.hh 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//	Add errors to a LorentzVector.  Used for direction errors 
//      BaBar native class
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Abi Soffer
//
// History
//      Mon Apr 27 00:28:13 PDT 1998, created from BbrVectorErr.hh, Abi Soffer
//
// Copyright Information:
//	Copyright (C) 1998
//
//------------------------------------------------------------------------
#ifndef BBRLORENTZVECTORERR_HH
#define BBRLORENTZVECTORERR_HH

#include <iosfwd>
#include <iosfwd>
class BbrVectorErr;


#include "BaBar/BaBar.hh"
#include "BbrGeom/BbrError.hh"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Geometry/Translation.h"
#include "CLHEP/Geometry/Transformation.h"

class BbrLorentzVectorErr : public HepLorentzVector {
public:
  enum MPolarCoordinateIndex {   
    Mom = 0, 
    Theta = 1, 
    Phi = 2,
    Mass = 3, 
    NUM_MPCOORDINATES = 4
  };
  enum EPolarCoordinateIndex {   
    Energy = 0, 
    ETheta = 1, 
    EPhi = 2,
    EMom = 3, 
    NUM_EPCOORDINATES = 4
  };
  // argumentless constructor:
  BbrLorentzVectorErr() : HepLorentzVector(), _covMatrix(NUM_COORDINATES) {}

  // auto casting constructor
  BbrLorentzVectorErr(const HepLorentzVector &p) : 
    HepLorentzVector(p), _covMatrix(NUM_COORDINATES)	{}

  BbrLorentzVectorErr(const HepLorentzVector &p, const BbrError& covMat) : 
    HepLorentzVector(p),
    _covMatrix(NUM_COORDINATES)				{ _covMatrix=covMat; }
  
  // Useful constructor for tracks:
  BbrLorentzVectorErr(const BbrVectorErr & p3, double mass);		      

  // copy constructor
  BbrLorentzVectorErr(const BbrLorentzVectorErr& v) : HepLorentzVector(v),
    _covMatrix(v.covMatrix())	{}

  // destructor MAY be needed later
  // virtual ~BbrLorentzVectorErr() {};

  // assignment operator:
  BbrLorentzVectorErr& operator=(const BbrLorentzVectorErr& v)
    {
      if (this != &v) {
	HepLorentzVector::operator=(v);
	_covMatrix = v.covMatrix();
      }
      return *this;
    }

  // mathematical modifiers:
  BbrLorentzVectorErr operator - () {
      HepLorentzVector t = *this;
      return BbrLorentzVectorErr( -t, _covMatrix);  // _covMatrix unaltered
  }

  BbrLorentzVectorErr& operator += (const BbrLorentzVectorErr& v){
      HepLorentzVector::operator+=(v);
      _covMatrix += v.covMatrix();
      return *this;
  }
  
  BbrLorentzVectorErr& operator -= (const BbrLorentzVectorErr& v){
      HepLorentzVector::operator-=(v);
      _covMatrix += v.covMatrix();
      return *this;
  }

  // can't implement this since there is no
  //  HepLorentzVector::transform(const HepTranslation):
  //
  //  BbrLorentzVectorErr& transform(const HepTranslation& trans){
  //    HepLorentzVector::transform(trans);
  //    return *this;
  //  }

  BbrLorentzVectorErr& transform(const HepRotation& rot);
  
  BbrLorentzVectorErr& transform(const HepLorentzRotation& rot);
  
  // can't implement this since there is no 
  // HepLorentzVector::transform(const HepTransformation):
  //
  //  BbrLorentzVectorErr& transform(const HepTransformation& transf){
  //    HepLorentzVector::transform(transf);
  //    _covMatrix = _covMatrix.similarity(transf.rot_mat());
  //    return *this;
  //  }

  double determineChisq(const HepLorentzVector& refVector) const;   
  // returns Chisquare
  // refVector refers to the same origin as the LorentzVector of this
  // ie refVector is not relative to this Vector

  BbrError covMRTPMatrix() const;
  // returns the covariance Matrix in spherical coordinate and mass
  // use   MPolarCoordinateIndex enum to get the components
  BbrError covETPRMatrix() const;
  // returns the covariance Matrix in spherical coordinate and mass
  // use   EPolarCoordinateIndex enum to get the components
  // note: it is different from the others because of the different EMC convention

  inline const BbrError& covMatrix() const    { return _covMatrix; }
  inline void setCovMatrix(const BbrError& v) { _covMatrix = v; }

//  void printOn(ostream& out=cout) const;

private:
  
  BbrError _covMatrix;
};

// globals:
BbrLorentzVectorErr operator + (const BbrLorentzVectorErr&, 
				const BbrLorentzVectorErr&);

BbrLorentzVectorErr operator - (const BbrLorentzVectorErr&, 
				const BbrLorentzVectorErr&);

std::ostream & operator<<(std::ostream & stream, const BbrLorentzVectorErr & verr);
std::istream & operator>>(std::istream & stream, BbrLorentzVectorErr & verr);

#endif






