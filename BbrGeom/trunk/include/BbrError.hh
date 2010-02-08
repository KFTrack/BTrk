//--------------------------------------------------------------------------
// File and Version Information:
//	$Id: BbrError.hh 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//      A wrapper for a covariance matrix.  A covariance matrix is
//	a symmetric n X n matrix.  Change in chisq from point
//	covariance matrix was determined is just
//
//	diff * covariance^-1 * diff_transpose
//
//	which is implemented in a similarity transform in HepSymMatrix
//	the method determineChisq carries this calculation out and requires
//	the result to be a scalar.
//      
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Forest Rouse                Jan 1996
//      Victoria Novotny            Aug 1996
//
// Copyright Information:
//	Copyright (C) 1996
//
//------------------------------------------------------------------------
#ifndef BBRERROR_HH
#define BBRERROR_HH

#include <iostream>
#include <math.h>

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/LorentzRotation.h"

class BbrError : public HepSymMatrix {

public:

  static const double chisqUndef;
  BbrError() : HepSymMatrix() {}

  BbrError(int n) : HepSymMatrix(n, 0)				{}

  // autocast copy constructor.  HepSymMatrix's promoted back
  // into BbrError matrices.

  BbrError( const HepSymMatrix &p ) : HepSymMatrix(p)		{}

  // new constructors for this class
  BbrError(const BbrError& v)					{*this = v;}

  BbrError& operator=(const BbrError& v)
    {
          if (this != &v) {
            HepSymMatrix::operator=(v);
      }
      return *this;
    }

  BbrError& operator=(const HepSymMatrix& v)
    {
          if (this != &v) {
            HepSymMatrix::operator=(v);
      }
      return *this;
    }

  // destructor MAY be needed later
  // virtual ~BbrError() {};

  //----------------------------------------------------------------------
  // determineChisq
  //	Compute v^T * V^(-1)*v - ie the chisq for this covariance
  //	matrix and the difference vector v.
  //----------------------------------------------------------------------

  double determineChisq(const HepVector& diff) const; 

  // Get right signature for all operations performed on BbrError matrices
  // that should (and will) result in BbrError matrices.  These include
  // similarity transforms, transpose, inverse, matrix multiplication,
  // addition, and subtraction.  HepSymMatrix's as a result of operations
  // are promoted back into BbrError matrices if we start out
  // with BbrError matrices in the first place. (See copy constructors)

  BbrError& operator *= (double t)
    { return (BbrError&) HepSymMatrix::operator*=(t); }

  BbrError& operator /= (double t)
    {return (BbrError&) HepSymMatrix::operator/=(t); }

  BbrError& operator += (const BbrError& m2)
    {
      HepSymMatrix::operator+=(m2);
      return *this;
    }

  BbrError& operator -= (const BbrError& m2)
    {
      HepSymMatrix::operator-=(m2);
      return *this;
    }

  BbrError operator - ()
    { BbrError temp(*this); 
      return temp; }
  // does nothing -- covariance Matrices have never negative entries on the 
  // main diagonal

  // Implement the signature for operators I also inherit
  BbrError& operator += (const HepSymMatrix& m2)
    {
      HepSymMatrix::operator+=(m2);
      return *this;
    }

  BbrError& operator -= (const HepSymMatrix& m2)
    {
      HepSymMatrix::operator-=(m2);
      return *this;
    }

  BbrError operator - () const
    { BbrError temp(*this); 
    return temp;}

  BbrError& operator += (const HepDiagMatrix& m2)
    {
      HepSymMatrix::operator+=(m2);
      return *this;
    }

  BbrError& operator -= (const HepDiagMatrix& m2)
    {
      HepSymMatrix::operator-=(m2);
      return *this;
    }

  BbrError similarity(const HepRotation& rot) const;
  BbrError similarity(const HepLorentzRotation& rot) const;
  // When feasible implement R * covMatrix * R_transpose (= R^-1)

  BbrError similarity(const BbrError& E);
  // implement E * covMatrix * E

  BbrError similarity(const HepMatrix& m1) const
    {
      BbrError mret(m1.num_row());
      mret.similarityWith(*this, m1);
      return mret;
    }

  BbrError& similarityWith(const BbrError& m, const HepMatrix& m1);

  // provide call ups to base classes similarity methods not implemented
  // here
  double similarity(const HepVector &v) const
       { return this->HepSymMatrix::similarity( v ); }
  HepSymMatrix similarity(const HepSymMatrix &m1) const
       { return this->HepSymMatrix::similarity( m1 ); }

private:
    
  friend BbrError operator*(double t, const BbrError& m1);

  friend BbrError operator*(const BbrError& m1, double t);

  friend BbrError operator/(double t, const BbrError& m1);

  friend BbrError operator/(const BbrError& m1, double t);

  friend BbrError operator+(const BbrError& m1, const BbrError& m2);
 
  friend BbrError operator-(const BbrError& m1, const BbrError& m2);
  
  friend std::ostream& operator<<(std::ostream& out, const BbrError& mat);
  friend std::istream& operator>>(std::istream& in, BbrError& mat);
  
};

#endif






