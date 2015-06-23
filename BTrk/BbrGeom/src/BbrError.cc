//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: BbrError.cc 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//      Class BbrError
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Forest Rouse            February 1996
//      Victoria Novotny        August   1996
//
// Copyright Information:
//      Copyright (C) 1996      U.C. Davis 
//
//------------------------------------------------------------------------
// File BbrError.cc
// Source file for class BbrError
//
// For advice, input, or any questions then please contact either
// Bob Jacobsen <Bob_Jacobsen@lbl.gov> or
// Forest Rouse <rouse@ucdhep.ucdavis.edu>
//
// =====================================================================
// Name           Change description
// Date
// Version
// =====================================================================
#include "BaBar/BaBar.hh"

static const char rscid[] = "$Id: BbrError.cc 491 2010-01-13 16:59:16Z stroili $";


#include <stdio.h>
#include <assert.h>
#include <ctype.h>


#include "BbrGeom/BbrError.hh"
#include "CLHEP/Matrix/Matrix.h"
using std::istream;
using std::ostream;

const double BbrError::chisqUndef = -1.;

BbrError BbrError::similarity(const HepRotation& rot) const
{
  HepMatrix mat(3,3); 
  mat(1,1)=rot.xx(); mat(1,2)=rot.xy(); mat(1,3)=rot.xz();
  mat(2,1)=rot.yx(); mat(2,2)=rot.yy(); mat(2,3)=rot.yz();
  mat(3,1)=rot.zx(); mat(3,2)=rot.zy(); mat(3,3)=rot.zz();

  HepSymMatrix w = similarity(mat);
  return w;
}


BbrError BbrError::similarity(const HepLorentzRotation& rot) const
{
  HepMatrix mat(4,4); 
  mat(1,1)=rot.xx(); mat(1,2)=rot.xy(); mat(1,3)=rot.xz(); mat(1,4)=rot.xt();
  mat(2,1)=rot.yx(); mat(2,2)=rot.yy(); mat(2,3)=rot.yz(); mat(2,4)=rot.yt();
  mat(3,1)=rot.zx(); mat(3,2)=rot.zy(); mat(3,3)=rot.zz(); mat(3,4)=rot.zt();
  mat(4,1)=rot.tx(); mat(4,2)=rot.ty(); mat(4,3)=rot.tz(); mat(4,4)=rot.tt();

  HepSymMatrix w = similarity(mat);
  return w;
}

BbrError BbrError::similarity(const BbrError& E)
{
  BbrError mret(HepSymMatrix::similarity(E));
  return mret;
}

BbrError& BbrError::similarityWith(const BbrError& mat,
                                           const HepMatrix& m1)
{
  assert(num_row() == m1.num_row());
  HepMatrix temp = m1*mat;
  register double tmp;

  for (int r = 0; r < num_row(); r++) {
    for (int c = 0; c <= r; c++) {
      tmp = 0.;
      for (int k = 0; k < m1.num_col(); k++) {
        tmp += temp[r][k]*m1[c][k];
      }
      (*this)[r][c] = tmp;
    }
  }
  return *this;
}




//----------------------------------------------------------------------
// determineChisq
//	Compute v^T * V^(-1)*v - ie the chisq for this covariance
//	matrix and the difference vector v.
//----------------------------------------------------------------------

double BbrError::determineChisq(const HepVector& diff) const 
{
  int ierr;
  HepMatrix dMat(diff.num_row(), 1);

  for (int i = 0; i < diff.num_row(); i++) dMat[i][0] = diff[i];
  
  double chisq = (inverse(ierr).similarityT(dMat))[0][0];

  if (ierr == 0) return chisq;
  else           return chisqUndef;
}

ostream& operator<<(ostream& out, const BbrError& mat)
{
    out << "Bbr Covariance Matrix:";
    out << (HepSymMatrix&) mat;
    return out;
}

istream& operator>>(istream& in, BbrError& mat) {
  // Peek at the next non-space character:
  char nextChar = ' ';
  while (isspace(nextChar)){
    nextChar = in.get();
  }
  in.putback(nextChar);

  if (EOF != nextChar){
    if (!isdigit(nextChar)){
      // Remove the "Bbr Covariance Matrix:" line:
      const int DUMMY_SIZE = 1000;
      char dummy[DUMMY_SIZE];
      in.getline(dummy, DUMMY_SIZE);
    }
    // Read in the matrix:
    for (int row = 1; row <= mat.num_row(); ++row){
      for (int col = 1; col <= mat.num_col(); ++col){
	in >> mat(row, col);
      }
    }
  }
  return in;
}


BbrError operator*(double t, const BbrError& m1)
  {
    BbrError mret = m1;
    mret *= t;
    return mret;
  }

BbrError operator*(const BbrError& m1, double t)
  {
    BbrError mret = m1;
    mret *= t;
    return mret;
  }

BbrError operator/(double t, const BbrError& m1)
  {
    BbrError mret = m1;
    mret /= t;
    return mret;
  }

BbrError operator/(const BbrError& m1, double t)
  {
    BbrError mret = m1;
    mret /= t;
    return mret;
  }

BbrError operator+(const BbrError& m1, const BbrError& m2)
  {
    BbrError mret = m1;
    mret += m2;
    return mret;
  }

BbrError operator-(const BbrError& m1, const BbrError& m2)
  {
    BbrError mret = m1;
    mret -= m2;
    return mret;
  }
  
 







