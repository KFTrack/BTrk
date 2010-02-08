//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: BbrLorentzVectorErr.cc 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//      Class BbrLorentzVectorErr
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Abi Soffer, Mon Apr 27 00:45:49 PDT 1998
//
// Copyright Information:
//      Copyright (C) 1996      U.C. Davis
//
//------------------------------------------------------------------------
// File BbrLorentzVectorErr.cc
// Source file for class BbrLorentzVectorErr
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
#include <float.h>
#include <iostream>
#include "difAlgebra/DifNumber.hh"
#include "difAlgebra/DifArray.hh"
#include "BbrGeom/BbrLorentzVectorErr.hh"
#include "BbrGeom/BbrVectorErr.hh"
using std::istream;
using std::ostream;



BbrLorentzVectorErr::BbrLorentzVectorErr(const BbrVectorErr & p3, 
					 double mass) :
  HepLorentzVector((const Hep3Vector &)p3, sqrt(p3.mag2() + mass * mass)),
  _covMatrix(NUM_COORDINATES)
{  
  // The 3-vector part of the error, initialized as a 0-matrix, then copied 
  // from p3:
  HepSymMatrix p4Err(4,0);
  int i, j;
  for (i = 1; i < 4; ++i){
    for (j = i; j < 4; ++j){  // start from j=i, since (i,j) = (j,i)
      p4Err(i,j) = p3.covMatrix()(i,j);
    }
  }

  // The energy part of the error:
  const double energy = t();

  if (energy != 0){
    p4Err(4,4) = 0.0;
    for (i = 0; i < 3; ++i){
      for (j = 0; j < 3; ++j){
	p4Err(4,4) += p3[i] * p3[j] *
	  p3.covMatrix()(i + 1, j + 1);
      }
    }
    p4Err(4,4) /= (energy * energy);
  }
  
  // The E-p correlated part of the error is
  // <(E - <E>) (Pi - <Pi>)> = dE/dPj Vjk dPi/dPk
  // Since dPi/dPk = delta_ik, we get dE/dPj Vji.
  for (i = 1; i < 4; ++i){
    for (j = 1; j < 4; ++j){
      p4Err(i,4) += p3(j-1) / energy * p3.covMatrix()(j,i);	
      // No need to switch indices for a HepSymMatrix, so no p4Err(4,i)
    }
  }

  setCovMatrix(p4Err);
}


 //-------------------------------------------------------------
BbrError BbrLorentzVectorErr::covMRTPMatrix() const{
  // protect against 0's
  double xv = x()==0 ?  FLT_MIN : x();
  double yv = y()==0 ?  FLT_MIN : y();
  double zv = z()==0 ?  FLT_MIN : z();
  double tv = t()==0 ?  FLT_MIN : t();
  DifNumber xDF(xv,X+1,NUM_MPCOORDINATES), yDF(yv,Y+1,NUM_MPCOORDINATES), zDF(zv,Z+1,NUM_MPCOORDINATES), tDF(tv,T+1,NUM_MPCOORDINATES);
  DifArray pars(NUM_MPCOORDINATES,NUM_COORDINATES);
  pars[Mom]   =  sqrt(xDF*xDF + yDF*yDF + zDF*zDF);
  pars[Mass]   = sqrt(fabs( tDF*tDF-(xDF*xDF + yDF*yDF + zDF*zDF)));
  pars[Phi]   = atan2(yDF,xDF);
  pars[Theta] = atan2(sqrt(xDF*xDF + yDF*yDF),zDF);
  BbrError result(covMatrix().similarity(pars.jacobian()));
  return result;
}

 //-------------------------------------------------------------
BbrError BbrLorentzVectorErr::covETPRMatrix() const{
  // protect against 0's
  double xv = x()==0 ?  FLT_MIN : x();
  double yv = y()==0 ?  FLT_MIN : y();
  double zv = z()==0 ?  FLT_MIN : z();
  double tv = t()==0 ?  FLT_MIN : t();
  DifNumber xDF(xv,X+1,NUM_EPCOORDINATES), yDF(yv,Y+1,NUM_EPCOORDINATES), zDF(zv,Z+1,NUM_EPCOORDINATES), tDF(tv,T+1,NUM_EPCOORDINATES);
  DifArray pars(NUM_EPCOORDINATES,NUM_COORDINATES);
  pars[EMom]   =  sqrt(xDF*xDF + yDF*yDF + zDF*zDF);
  pars[Energy]   =  tDF;
  pars[EPhi]   = atan2(yDF,xDF);
  pars[ETheta] = atan2(sqrt(xDF*xDF + yDF*yDF),zDF);
  BbrError result(covMatrix().similarity(pars.jacobian()));
  return result;
}


//-------------------------------------------------------------
double 
BbrLorentzVectorErr::determineChisq(const HepLorentzVector& refVector) const
{
   HepVector temp(NUM_COORDINATES, 0);
   temp[0] = refVector.x()-this->x();
   temp[1] = refVector.y()-this->y();
   temp[2] = refVector.z()-this->z();
   temp[3] = refVector.t()-this->t();

   return _covMatrix.determineChisq(temp);
}


//-------------------------------------------------------------
BbrLorentzVectorErr& 
BbrLorentzVectorErr::transform(const HepLorentzRotation& rot) {
  HepLorentzVector::transform(rot);
  _covMatrix = _covMatrix.similarity(rot);
  return *this;
}

//-------------------------------------------------------------
BbrLorentzVectorErr& 
BbrLorentzVectorErr::transform(const HepRotation& rot){
  HepLorentzVector::transform(rot);
  HepMatrix tempRot(NUM_COORDINATES, NUM_COORDINATES);

  // Fill a 4x4 matrix from the 3x3 HepRotation. Note that they use different
  // indexing schemes (!?@#$&^*&#$@#):

  int row;
  int col;
  for (row = 0; row < 3; ++row){   // 3 is the size of HepRotation
    for (col = 0; col < 3; ++col){ // (which provides no enum)
      tempRot(row+1, col+1) = rot(row, col);
    }
  }
  
  // fill the 4th row:
  tempRot(4,4) = 1.0;
  for (col = 1; col < 4; ++col){
    tempRot(4, col) = 0.0;
  }

  // fill the 4th column:
  for (row = 1; row < 4; ++row){
    tempRot(row, 4) = 0.0;
  }

  _covMatrix = _covMatrix.similarity(tempRot);
  return *this;
}



//-------------------------------------------------------------
BbrLorentzVectorErr 
operator + (const BbrLorentzVectorErr& v, const BbrLorentzVectorErr& w){
  BbrLorentzVectorErr ve(HepLorentzVector(v.x()+w.x(),v.y()+w.y(),v.z()+w.z(),
					  v.t()+w.t()),
			 (v.covMatrix()+w.covMatrix()));
  return ve;
}


//-------------------------------------------------------------
BbrLorentzVectorErr 
operator - (const BbrLorentzVectorErr& v, const BbrLorentzVectorErr& w){
  BbrLorentzVectorErr ve(HepLorentzVector(v.x()-w.x(),v.y()-w.y(),v.z()-w.z(),
					  v.t()-w.t()), 
			 (v.covMatrix()+w.covMatrix()));
  return ve;
}



//-------------------------------------------------------------
ostream & operator<<(ostream & stream, const BbrLorentzVectorErr & verr) {
  stream << (const HepLorentzVector&)verr
	 << ", " << verr.covMatrix();
  
  return stream;
}

istream & operator>>(istream & stream, BbrLorentzVectorErr & verr) {
  BbrError mat(verr.SIZE);
  stream >> (HepLorentzVector&)verr >> mat;
  verr.setCovMatrix(mat);
  
  return stream;
}
