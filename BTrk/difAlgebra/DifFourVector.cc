//--------------------------------------------------------------------------
// File and Version Information: 
// 	$Id: DifFourVector.cc 501 2010-01-14 12:46:50Z stroili $
//
// Description:
//	Class Implementation for |DifFourVector| 
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

#include "BTrk/difAlgebra/DifFourVector.hh"
using std::endl;
using std::ostream;
using namespace CLHEP;


DifFourVector::DifFourVector()
  :E(0.0),P(0.0,0.0,0.0)
{} 

DifFourVector::DifFourVector
(const DifNumber& m,const DifVector& p)
  :E(sqrt(m*m+p*p)),P(p)
{}
DifFourVector::DifFourVector
(const double& m,const DifVector& p)
  :E(sqrt(m*m+p*p)),P(p)
{}

DifFourVector::DifFourVector(const DifFourVector& v)
  :E(v.E),P(v.P)
{}

HepSymMatrix DifFourVector::errorMatrix(const HepSymMatrix& e)const {
  HepSymMatrix temp(4);
  temp(1,1)=correlation(E,E,e);  
  temp(1,2)=correlation(E,P.x,e);
  temp(1,3)=correlation(E,P.y,e);
  temp(1,4)=correlation(E,P.z,e);
  temp(2,2)=correlation(P.x,P.x,e);
  temp(2,3)=correlation(P.x,P.y,e);
  temp(2,4)=correlation(P.x,P.z,e);
  temp(3,3)=correlation(P.y,P.y,e);
  temp(3,4)=correlation(P.y,P.z,e);
  temp(4,4)=correlation(P.z,P.z,e);
  return temp;
}

HepMatrix DifFourVector::jacobian()const{
  int npar=E.nPar();
  HepMatrix temp(4,npar);
  for(int i=1; i<=npar; i++){
    temp(1,i)=E.derivative(i);
    temp(2,i)=P.x.derivative(i);
    temp(3,i)=P.y.derivative(i);
    temp(4,i)=P.z.derivative(i);
  } // (int i=1; i<=npar; i++)
  return temp;
}

void
DifFourVector::boostTo(const DifFourVector& pTo)
{
  const DifVector xHat(1,0,0);
  const DifVector yHat(0,1,0);
  const DifVector zHat(0,1,0);
  DifVector z=pTo.direction();
  DifVector y=zHat-z*(zHat*z);
  if(y.length()<0.01) y=xHat-z*(xHat*z);
  y.normalize();
  DifVector x(cross(y,z));
  
  DifNumber px=P*x;
  DifNumber py=P*y;
  DifNumber pz=P*z;

  DifNumber gamma=pTo.E/pTo.mass();
  DifNumber beta=pTo.pMag()/pTo.E;

  DifNumber pzP=gamma*(pz-beta*E);
  DifNumber eP=gamma*(E-beta*pz);

  E=eP;
  P=px*x+py*y+pzP*z;

  return;

}

void
DifFourVector::boostToMe(std::vector<DifFourVector*>& list)const{
  const DifVector xHat(1,0,0);
  const DifVector yHat(0,1,0);
  const DifVector zHat(0,0,1);
  DifVector z=P;
  z.normalize();
  DifVector y(zHat-z*(zHat*z));
  if(y.lengthSq()<0.0001) y=xHat-z*(xHat*z);
  y.normalize();
  DifVector x(cross(y,z));
  
  DifNumber gamma=E/mass();
  DifNumber beta=pMag()/E;
  
  for(size_t i=0;i<list.size();i++) {
    DifFourVector& p4i=*list[i];
    DifNumber px=p4i.P*x;
    DifNumber py=p4i.P*y;
    DifNumber pz=p4i.P*z;
    DifNumber e=p4i.E;

    DifNumber pzP=gamma*(pz-beta*e);
    DifNumber eP=gamma*(e-beta*pz);
    
    p4i.E=eP;
    p4i.P=px*x+py*y+pzP*z;

  }
    
}


void
DifFourVector::boostFromMe(std::vector<DifFourVector*>& list)const{
  const DifVector xHat(1,0,0);
  const DifVector yHat(0,1,0);
  const DifVector zHat(0,0,1);
  DifVector z=P;
  z.normalize();
  DifVector y(zHat-z*(zHat*z));
  if(y.lengthSq()<0.0001) y=xHat-z*(xHat*z);
  y.normalize();
  DifVector x(cross(y,z));
  
  DifNumber gamma=E/mass();
  DifNumber beta=pMag()/E;
  
  for(size_t i=0;i<list.size();i++) {
    DifFourVector& p4i=*list[i];
    DifNumber px=p4i.P*x;
    DifNumber py=p4i.P*y;
    DifNumber pz=p4i.P*z;
    DifNumber e=p4i.E;

    DifNumber pzP=gamma*(pz+beta*e);
    DifNumber eP=gamma*(e+beta*pz);
    
    p4i.E=eP;
    p4i.P=px*x+py*y+pzP*z;

  }
    
}

void
DifFourVector::boostFrom(const DifFourVector& pFrom)
{
 const DifVector xHat(1,0,0);
  const DifVector yHat(0,1,0);
  const DifVector zHat(0,1,0);
  DifVector z=pFrom.direction();
  DifVector y=zHat-z*(zHat*z);
  if(y.length()<0.01) y=xHat-z*(xHat*z);
  y.normalize();
  DifVector x(cross(y,z));
  
  DifNumber px=P*x;
  DifNumber py=P*y;
  DifNumber pz=P*z;

  DifNumber gamma=pFrom.E/pFrom.mass();
  DifNumber beta=pFrom.pMag()/pFrom.E;

  DifNumber pzP=gamma*(pz+beta*E);
  DifNumber eP=gamma*(E+beta*pz);

  E=eP;
  P=px*x+py*y+pzP*z;
}


void DifFourVector::print(ostream& o)const {
  o << "E:" << endl << E;
  o << "P:" << endl << P;
}

