//--------------------------------------------------------------------------
// File and Version Information: 
// 	$Id: DifVector.cc 501 2010-01-14 12:46:50Z stroili $
//
// Description:
//	Class Implementation for |DifVector| 
//      What do i do ?
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
#include "BTrk/BaBar/BaBar.hh"

#include "BTrk/difAlgebra/DifVector.hh"
#include "BTrk/difAlgebra/DifRotation.hh"
#include "CLHEP/Vector/ThreeVector.h"
using std::ostream;
using namespace CLHEP;

extern const DifVector xhat(1,0,0);
extern const DifVector yhat(0,1,0);
extern const DifVector zhat(0,0,1);
extern const DifVector nullVec(0,0,0);

DifVector::DifVector()
  :x(0.0),y(0.0),z(0.0)  
{} 

DifVector::DifVector
(const DifNumber& X,const DifNumber& Y,const DifNumber& Z)
  :x(X),y(Y),z(Z)
{}

DifVector::DifVector
(double X, double Y,double Z)
  :x(X),y(Y),z(Z)
{}

DifVector::DifVector(const Hep3Vector& v)
  :x(v.x()),y(v.y()),z(v.z())
{}

DifVector::DifVector(const DifVector& v)
  :x(v.x),y(v.y),z(v.z)
{}

HepSymMatrix DifVector::errorMatrix(const HepSymMatrix& e)const {
  HepSymMatrix temp(3);
  temp(1,1)=correlation(x,x,e);  
  temp(1,2)=correlation(x,y,e);
  temp(1,3)=correlation(x,z,e);
  temp(2,2)=correlation(y,y,e);
  temp(2,3)=correlation(y,z,e);
  temp(3,3)=correlation(z,z,e);
  return temp;
}

HepMatrix DifVector::jacobian()const{
  int npar=x.nPar();
  HepMatrix temp(3,npar);
  for(int i=1; i<=npar; i++){
    temp(1,i)=x.derivative(i);
    temp(2,i)=y.derivative(i);
    temp(3,i)=z.derivative(i);
  } // (int i=1; i<=npar; i++)
  return temp;
}

DifVector& DifVector::rotate(const DifRotation& r) {
  r.rotate(*this); return *this;
}

DifVector& DifVector::rotate(const DifNumber& a,const DifNumber& b,const DifNumber& g) {
  DifRotation(a,b,g).rotate(*this);
  return *this;
}

void DifVector::print(ostream& o)const {
  o << "x:" << x;
  o << "y:" << y;
  o << "z:" << z;
}
