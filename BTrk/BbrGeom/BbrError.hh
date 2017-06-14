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
//	which is implemented in a similarity transform in CLHEP::HepSymMatrix
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

#include <math.h>
#include <iostream>

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/Rotation.h"

class BbrError : public CLHEP::HepSymMatrix {
   public:
    static const double chisqUndef;
    BbrError() : CLHEP::HepSymMatrix() {}

    BbrError(int n) : CLHEP::HepSymMatrix(n, 0) {}

    // autocast copy constructor.  CLHEP::HepSymMatrix's promoted back
    // into BbrError matrices.

    BbrError(const CLHEP::HepSymMatrix& p) : CLHEP::HepSymMatrix(p) {}

    // new constructors for this class
    BbrError(const BbrError& v) { *this = v; }

    BbrError& operator=(const BbrError& v) {
        if (this != &v) {
            CLHEP::HepSymMatrix::operator=(v);
        }
        return *this;
    }

    BbrError& operator=(const CLHEP::HepSymMatrix& v) {
        if (this != &v) {
            CLHEP::HepSymMatrix::operator=(v);
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

    double determineChisq(const CLHEP::HepVector& diff) const;

    // Get right signature for all operations performed on BbrError matrices
    // that should (and will) result in BbrError matrices.  These include
    // similarity transforms, transpose, inverse, matrix multiplication,
    // addition, and subtraction.  CLHEP::HepSymMatrix's as a result of operations
    // are promoted back into BbrError matrices if we start out
    // with BbrError matrices in the first place. (See copy constructors)

    BbrError& operator*=(double t) { return (BbrError&)CLHEP::HepSymMatrix::operator*=(t); }

    BbrError& operator/=(double t) { return (BbrError&)CLHEP::HepSymMatrix::operator/=(t); }

    BbrError& operator+=(const BbrError& m2) {
        CLHEP::HepSymMatrix::operator+=(m2);
        return *this;
    }

    BbrError& operator-=(const BbrError& m2) {
        CLHEP::HepSymMatrix::operator-=(m2);
        return *this;
    }

    BbrError operator-() {
        BbrError temp(*this);
        return temp;
    }
    // does nothing -- covariance Matrices have never negative entries on the
    // main diagonal

    // Implement the signature for operators I also inherit
    BbrError& operator+=(const CLHEP::HepSymMatrix& m2) {
        CLHEP::HepSymMatrix::operator+=(m2);
        return *this;
    }

    BbrError& operator-=(const CLHEP::HepSymMatrix& m2) {
        CLHEP::HepSymMatrix::operator-=(m2);
        return *this;
    }

    BbrError operator-() const {
        BbrError temp(*this);
        return temp;
    }

    BbrError& operator+=(const CLHEP::HepDiagMatrix& m2) {
        CLHEP::HepSymMatrix::operator+=(m2);
        return *this;
    }

    BbrError& operator-=(const CLHEP::HepDiagMatrix& m2) {
        CLHEP::HepSymMatrix::operator-=(m2);
        return *this;
    }

    BbrError similarity(const CLHEP::HepRotation& rot) const;
    BbrError similarity(const CLHEP::HepLorentzRotation& rot) const;
    // When feasible implement R * covMatrix * R_transpose (= R^-1)

    BbrError similarity(const BbrError& E);
    // implement E * covMatrix * E

    BbrError similarity(const CLHEP::HepMatrix& m1) const {
        BbrError mret(m1.num_row());
        mret.similarityWith(*this, m1);
        return mret;
    }

    BbrError& similarityWith(const BbrError& m, const CLHEP::HepMatrix& m1);

    // provide call ups to base classes similarity methods not implemented
    // here
    double similarity(const CLHEP::HepVector& v) const {
        return this->CLHEP::HepSymMatrix::similarity(v);
    }
    CLHEP::HepSymMatrix similarity(const CLHEP::HepSymMatrix& m1) const {
        return this->CLHEP::HepSymMatrix::similarity(m1);
    }

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
