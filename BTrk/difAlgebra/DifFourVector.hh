//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DifFourVector.hh 501 2010-01-14 12:46:50Z stroili $
//
// Description:
//	Class Header for |DifVector|
//      A 4-vector based on differential numbers
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	A. Snyder
//
// Copyright Information:
//	Copyright (C) 2002	SLAC
//
//------------------------------------------------------------------------

#ifndef DifFourVector_HH
#define DifFourVector_HH

#include <vector>

#include "BTrk/difAlgebra/DifNumber.hh"
#include "BTrk/difAlgebra/DifVector.hh"

#include <iosfwd>
class DifRotation;
#include "CLHEP/Vector/ThreeVector.h"

class DifFourVector {

public:

  //constructors
  DifFourVector();			// null - default
  DifFourVector			// construct from components
  (const DifNumber& mass,const DifVector& p);
  DifFourVector			// construct from components
  (const double& mass,const DifVector& p);
  
  DifFourVector(const DifFourVector& v); // copy

  //destructor
  ~DifFourVector() {}		// destroy

  //given error on parameters calculate error on vector
  CLHEP::HepSymMatrix errorMatrix		// caclulate error matrix
  (const CLHEP::HepSymMatrix& e)const;
  
  //give jacobian transformation matrix wrt independent parameters;
  CLHEP::HepMatrix jacobian()const;

  //operators 
  inline DifFourVector& operator+=(const DifFourVector& v);
  inline DifFourVector& operator-=(const DifFourVector& v);
  inline DifFourVector& operator=(const DifFourVector & v);
  inline DifFourVector operator-()const;

  inline friend DifFourVector operator+	// vector sum
  (const DifFourVector& a,const DifFourVector& b);
  inline friend DifFourVector operator-	// vector difference
  (const DifFourVector& a,const DifFourVector& b); 
  inline friend DifNumber operator*	// scalar product
  (const DifFourVector& a,const DifFourVector& b);



  //access 
  inline int nPar()const;	// return number of params
  inline DifVector direction()const;
  inline DifNumber pMag()const {return P.length();}
  inline DifNumber massSq()const {return E*E-P*P;}
  inline DifNumber mass()const {
    DifNumber temp=massSq();
    if(temp>=0) return sqrt(temp);
    return -sqrt(-massSq());
  }


  //i/o
  void print(std::ostream& o)const;	// print out


  //manipulations

  //misc
  inline DifFourVector& zeroDerivatives(); // zero derivatives
  

  //boost


  void boostTo(const DifFourVector&);
  void boostFrom(const DifFourVector&);

  void boostToMe
  (std::vector<DifFourVector*>& listToBoost)const;
  void boostFromMe
  (std::vector<DifFourVector*>& listToBoost)const;
    

  //algebra


  //data members - public .. yes, folks that's intentional!
public:

  // energy-momentum components of a 4-vector

  DifNumber E;			// energy-like component
  DifVector P;			// momentum-like compoent


};

//io 
inline std::ostream& operator<<(std::ostream& o,const DifFourVector& n) {
  n.print(o);
  return o;
}

#include "BTrk/difAlgebra/DifFourVector.icc"

#endif









