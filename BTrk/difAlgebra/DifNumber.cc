//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DifNumber.cc 501 2010-01-14 12:46:50Z stroili $
//
// Description:
//	Class Implementation for |DifNumber| 
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

#include "BTrk/difAlgebra/DifNumber.hh"
#include "BTrk/difAlgebra/DifIndepPar.hh"
using std::endl;
using std::ostream;
using namespace CLHEP;

extern const DifNumber zero(0.0);
extern const DifNumber one(1.0);


double DifNumber::error(const HepSymMatrix& e)const {
  return sqrt(correlation(*this,e));
}

double DifNumber::error()const{
  if(indepPar()==0) {return 0.0;}
  return error(indepPar()->covariance());
}

void DifNumber::fetchDerivatives(HepVector& v)const {
  assert(v.num_row()==nPar());
  for(int i=1; i<=nPar(); i++) {v(i)=derivative(i);}
}

HepVector DifNumber::derivatives()const{
    HepVector temp(nPar());
    fetchDerivatives(temp);
    return temp;
}

double 
DifNumber::correlation(const DifNumber& b,const HepSymMatrix &e) const {
  assert(e.num_col()==nPar());
  assert(e.num_row()==b.nPar());
  double error = 0.;
  for(int i=1; i<=nPar(); i++) {
    for(int j=1; j<=b.nPar(); j++) {
      error+=derivative(i)*e(i,j)*b.derivative(j);
    }
  }
  return error;
}

double
DifNumber::correlation(const DifNumber& b)const {
  if(indepPar()==0) return 0.0;
  if(b.indepPar()!=indepPar()) return 0.0;
  return correlation(b,indepPar()->covariance());
}

double correlation(const DifNumber& a,const DifNumber& b) // correlation from default matrix
{ 
  return (a.indepPar()==0||b.indepPar()==0||a.indepPar()!=b.indepPar())?0:a.correlation(b,a.indepPar()->covariance()); 
}

// FIXME: This function should be inlined, but that would require checking out additional packages...
double correlation(const DifNumber& a,const DifNumber& b,const HepSymMatrix& e)	// correlation for specified error
{ return a.correlation(b,e); }

void DifNumber::print(ostream& o)const {
  o << "number:" << number() << endl;
  o << "npar:" << nPar() << endl;
  for(int i=1; i<=nPar(); i++) {
    o << "derivative(" << i << "):" << derivative(i) << endl;
  }
}

 
extern DifNumber solveQuad
(const DifNumber& a,		// quadratic term
 const DifNumber& b,		// linear term
 const DifNumber& c,		// const term
 int pref,			// solution preference
 Code& code)			// error code
{
  DifNumber descr=b*b-4.0*a*c;
  if(descr<0.0) {		// solution not real
    code.setFail(1341);
    return DifNumber(0.0);
  }
  if(a.number()==0.0){
    if(b.number()==0.0) {
      code.setFail(1342);
      return DifNumber(0.0);
    }
    code.setSuccess(40);
    return -c/b+a*c/pow(b,3);
  }
  code.setSuccess(40);
  descr=sqrt(descr);
  DifNumber s=-b;

  if(pref==+1) {		// positive solution
    s+=descr;
  }else if(pref==-1){		// negative solution
    s-=descr;
  }else if(pref==0) {		// smallest solution
    if(s>0.0) {s-=descr;}else {s+=descr;}
  }else {			// illegal prefrence
    code.setFail(1343);
    return DifNumber(0.0);
  }
  s/=2.0*a;
  return s;
}

