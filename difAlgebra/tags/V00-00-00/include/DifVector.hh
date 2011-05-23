//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DifVector.hh 501 2010-01-14 12:46:50Z stroili $
//
// Description:
//	Class Header for |DifVector|
//      A 3-vector based on differential numbers
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	A. Snyder
//
// Copyright Information:
//	Copyright (C) 1996	SLAC
//
//------------------------------------------------------------------------

#ifndef DifVector_HH
#define DifVector_HH

#include "difAlgebra/DifNumber.hh"

#include <iosfwd>
class DifRotation;
#include "CLHEP/Vector/ThreeVector.h"
#include <assert.h>

class DifVector {

public:

  //constructors
  DifVector();			// null - default
  DifVector			// construct from components
  (const DifNumber& x,const DifNumber& y,const DifNumber& z);
  DifVector
  (double x,double y,double z);
  DifVector			// construct from |Hep3Vector|
  (const Hep3Vector& v);
  DifVector(const DifVector& v); // copy

  //destructor
  ~DifVector() {}		// destroy

  //given error on parameters calculate error on vector
  HepSymMatrix errorMatrix		// caclulate error matrix
  (const HepSymMatrix& e)const;
  
  //give jacobian transformation matrix wrt independent parameters;
  HepMatrix jacobian()const;

  //operators 
  inline DifVector& operator+=(const DifVector& v);
  inline DifVector& operator-=(const DifVector& v);
  inline DifVector& operator*=(const DifNumber& n);
  inline DifVector& operator*=(const double& n);
  DifVector& operator/=(const DifNumber& n);
  DifVector& operator/=(const double& n);
  inline DifVector& operator=(const DifVector & v);
  inline DifVector operator-()const;

  inline friend DifVector operator+	// vector sum
  (const DifVector& a,const DifVector& b);
  inline friend DifVector operator-	// vector difference
  (const DifVector& a,const DifVector& b); 
  inline friend DifNumber operator*	// scalar product
  (const DifVector& a,const DifVector& b);
  inline friend DifVector operator*	// vector*scalar
  (const DifVector& v,const DifNumber& n);
  inline friend DifVector operator*	// scalar*vector
  (const DifNumber& n,const DifVector& v);
  inline friend DifVector operator/	// vector/scalar
  (const DifVector& v,const DifNumber& n);
  inline friend DifVector operator*	// vector*scalar
  (const DifVector& v,const double& n);
  inline friend DifVector operator*	// scalar*vector
  (const double& n,const DifVector& v);
  inline friend DifVector operator/	// vector/scalar
  (const DifVector& v,const double& n);


  //other operations 
  inline friend DifVector cross(const DifVector& a,const DifVector& b);

  //access 
  inline int nPar()const;	// return number of params

  //i/o
  void print(std::ostream& o)const;	// print out


  //manipulations

  //misc
  inline DifVector& flipsign();	// flip sign of all components
  inline DifVector& normalize();	// norm to unit vector
  inline DifVector& zeroDerivatives(); // zero derivatives
  

  //rotations
  inline DifVector& rotate(const DifVector& axis,const DifNumber& angle); // rotate around axis
  inline DifVector& rotate(const DifVector& axis, const DifNumber& cosine, const DifNumber& sine);// rotate with cos and sine
  DifVector& rotate(const DifNumber& alpha,const DifNumber& beta,const DifNumber& gamma); // euler angles
  DifVector& rotate(const DifRotation& r); // rotatation matrix
  inline DifVector& rotateX(const DifNumber& angle);	// around x
  inline DifVector& rotateX(const DifNumber& cosine,const DifNumber& sine);
  inline DifVector& rotateY(const DifNumber& angle);	// around y
  inline DifVector& rotateY(const DifNumber& cosine,const DifNumber& sine);
  inline DifVector& rotateZ(const DifNumber& angle);	// around z
  inline DifVector& rotateZ(const DifNumber& cosine,const DifNumber& sine);

  inline DifVector& rotate(const DifVector& axis,const double& angle); // rotate around axis
  inline DifVector& rotate(const DifVector& axis, const double& cosine, const double& sine);// rotate with cos and sine
  DifVector& rotate(const double& alpha,const double& beta,const double& gamma); // euler angles
  inline DifVector& rotateX(const double& angle);	// around x
  inline DifVector& rotateX(const double& cosine,const double& sine);
  inline DifVector& rotateY(const double& angle);	// around y
  inline DifVector& rotateY(const double& cosine,const double& sine);
  inline DifVector& rotateZ(const double& angle);	// around z
  inline DifVector& rotateZ(const double& cosine,const double& sine);

  //algebra
  inline  DifVector transverse	// part tranverse to |v|
  (const DifVector& v)const;
  inline  DifNumber dot		// scalr product
  (const DifVector& v)const; 
  inline DifNumber length()const; // length of vector
  inline DifNumber lengthSq()const; // length squared
  inline DifVector unit()const;	// direction
  inline DifNumber perp()const; // perp comp
  inline DifNumber perpSq()const; // perp squared


  //polar corrdinates
  DifNumber r()const;		// length by any other name
  DifNumber phi()const;		// azimutal angle
  DifNumber theta()const;	// polar angle
  DifNumber cosTheta()const;	// cosine of polar angle

  //data members - public .. yes, folks that's intentional!
public:

  //x,y,z components of 3-vector
  DifNumber x;
  DifNumber y;
  DifNumber z;


};

//io 
inline std::ostream& operator<<(std::ostream& o,const DifVector& n) {
  n.print(o);
  return o;
}

#include "difAlgebra/DifVector.icc"

#endif









