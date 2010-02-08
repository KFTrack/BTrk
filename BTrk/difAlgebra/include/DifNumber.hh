//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DifNumber.hh 501 2010-01-14 12:46:50Z stroili $
//
// Description:
//	Class Header for |DifNumber|
//      class to represent a number with derivatives
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

#ifndef DifNumber_HH
#define DifNumber_HH


#define MATRIX_BOUND_CHECK
#include "BaBar/BaBar.hh"
#include <math.h>
#include <assert.h>
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "difAlgebra/Code.hh"


#include <iosfwd>
#include <assert.h>

class DifIndepPar;
 
class DifNumber {

private:
  enum{_nmax=100};

public:

  //constructors and destructor
  inline
  DifNumber();			// zero with derivatives=0 
  inline
  explicit 			// disallow implicit promotion to DifNumber
  DifNumber			// s w/o derivatives 
  (double s,int npar=0);	// or s with derivatives=0
  inline
  DifNumber			// s with derivatives d
  (double s,const HepVector& d); 
  inline
  DifNumber(double s,int i,double di,int npar);         // s with i'th derivative set
  inline
  DifNumber(double s,int i,int npar);                   // s with ith derivative=1.0
  inline
  DifNumber(double s,const DifIndepPar* indeppar);      // s with derivatives 0 wrt |indeppar|
  inline
  DifNumber(double s,int i,const DifIndepPar* indepar); // s with i'th derivative wrt |indpar|
  
private:
  inline
  DifNumber(double newval,const DifNumber& old,double factor=0.0);// new value, derivs*factor

public:
  inline
  DifNumber(const DifNumber& s); // copy
  virtual ~DifNumber() {}

  //set
  inline void setNumber(double i)                  { _number=i;}
  inline void setNPar(int i)                       { _npar=i;}
  inline void setDerivatives(const HepVector& d)   { copyDerivs(d); }
  inline void setDerivative(int i,double value)    { _derivatives[i-1]=value;}
  inline void zeroDerivatives()                    { for(int i=1; i<=nPar(); i++) setDerivative(i,0.0); }
  inline void setIndepPar(const DifIndepPar* par)  { _indepPar=par;}
  inline void extendPar(int nnew);

  //access info
  inline double number()const {return _number;}
  inline double& gimeNumber() {return _number;}
  HepVector derivatives()const;
  inline double derivative(int i)const {return _derivatives[i-1];}
  inline int nMax() {return _nmax;}
  inline int nPar()const {return _npar;}
  inline const DifIndepPar* indepPar()const {return _indepPar;}
  double error(const HepSymMatrix& e)const; // error from specified error
  double error()const;		// error from default error
  inline void tickle()const {return;} // tickle a number

  //io
  void print(std::ostream& o)const;	// print out

  //operators and such

  //internals
  inline DifNumber& operator+=(const DifNumber& a);
  inline DifNumber& operator-=(const DifNumber& a);
  inline DifNumber& operator*=(const DifNumber& a);
  inline DifNumber& operator/=(const DifNumber& a);
  inline DifNumber& operator=(const DifNumber& a);
  inline DifNumber operator-()const 
  { return DifNumber(*this).flipsign(); } // unary minus

  inline DifNumber& operator+=(const double& a);
  inline DifNumber& operator-=(const double& a);
  inline DifNumber& operator*=(const double& a);
  inline DifNumber& operator/=(const double& a);
  inline DifNumber& operator=(const double& a);

  
  //manupulate this guy - reverse polish style algebra
  inline DifNumber& inverse();	// invert
  inline DifNumber& flipsign();	// flip the sign
  inline DifNumber& absolute()	/* take absolute value */      { if (_number<0) flipsign(); return *this; }
  inline void cosAndSin(DifNumber& c,DifNumber& s)const;	// calc sine and cosine of this
  inline DifNumber& squareRoot();		// take square root of this
  inline DifNumber& mod(double lo,double hi);	// mod to specified range
  inline DifNumber& arcTangent(const DifNumber& x); // arctangent |*this/x| = atan2(*this,x)
  inline DifNumber& power(double p); // raise to a power
  inline DifNumber& power(const DifNumber& p); // raise to a power

  inline double sign()const	/* return sign of number */  { return _number>=0 ? 1 : -1; }
  double correlation(const DifNumber& b,const HepSymMatrix& e) const;	// correlation for specified error
  double correlation(const DifNumber& b) const; // correlation from default matrix

  //sum into functions
  void sumMatrix(HepMatrix& m)const;

  //fetch functions
  void fetchNumber(double& n)const  {n=number();}
  void fetchDerivatives(HepVector& v)const;
  
  //friends
  inline friend DifNumber operator+(const DifNumber& a,const DifNumber& b) { return DifNumber(a)+=b; }
  inline friend DifNumber operator-(const DifNumber& a,const DifNumber& b) { return DifNumber(a)-=b; }
  inline friend DifNumber operator*(const DifNumber& a,const DifNumber& b) { return DifNumber(a)*=b; }
  inline friend DifNumber operator/(const DifNumber& a,const DifNumber& b) { return DifNumber(a)/=b; }

  inline friend DifNumber operator+(const DifNumber& a,const double& b) { return DifNumber(a)+=b; }
  inline friend DifNumber operator-(const DifNumber& a,const double& b) { return DifNumber(a)-=b; }
  inline friend DifNumber operator*(const DifNumber& a,const double& b) { return DifNumber(a)*=b; }
  inline friend DifNumber operator/(const DifNumber& a,const double& b) { return DifNumber(a)/=b; }

  inline friend DifNumber operator+(const double& a,const DifNumber& b) { return DifNumber(a)+=b; }
  inline friend DifNumber operator-(const double& a,const DifNumber& b) { return DifNumber(a)-=b; }
  inline friend DifNumber operator*(const double& a,const DifNumber& b) { return DifNumber(a)*=b; }
  inline friend DifNumber operator/(const double& a,const DifNumber& b) { return DifNumber(a)/=b; }


  inline friend bool operator>(const DifNumber& a,const DifNumber& b)  { return a.number()>b.number(); }
  inline friend bool operator>(const DifNumber& a,const double& b)  { return a.number()>b; }
  inline friend bool operator>(const double& a,const DifNumber& b)  { return a>b.number(); }

  inline friend bool operator<(const DifNumber& a,const DifNumber& b)  { return b>a; }
  inline friend bool operator<(const DifNumber& a,const double& b)  { return b>a; }
  inline friend bool operator<(const double& a,const DifNumber& b)  { return b>a; }


  inline friend bool operator>=(const DifNumber& a,const DifNumber& b) { return a.number()>=b.number(); }
  inline friend bool operator>=(const DifNumber& a,const double& b) { return a.number()>=b; }
  inline friend bool operator>=(const double& a,const DifNumber& b) { return a>=b.number(); }

  inline friend bool operator<=(const DifNumber& a,const DifNumber& b) { return b>=a; }
  inline friend bool operator<=(const DifNumber& a,const double& b) { return b>=a; }
  inline friend bool operator<=(const double& a,const DifNumber& b) { return b>=a; }

  inline friend bool operator==(const DifNumber& a,const DifNumber& b);
  inline friend bool operator==(const DifNumber& a,const double & b) {return false;}
  inline friend bool operator==(const double& a,const DifNumber& b) {return false;}

  inline friend bool operator!=(const DifNumber& a,const DifNumber& b) { return !(a==b); }
  inline friend bool operator!=(const DifNumber& a,const double& b) { return true; }
  inline friend bool operator!=(const double& a,const DifNumber& b) { return true; }
  
  inline friend DifNumber sin(const DifNumber& a)                      
  { return DifNumber(sin(a.number()),a,cos(a.number())); }
  inline friend DifNumber cos(const DifNumber& a)                     
  { return DifNumber(cos(a.number()),a,-sin(a.number())); }
  inline friend DifNumber tan(const DifNumber &a)                      
  { double t=tan(a.number()); return DifNumber(t,a,1.0+t*t); }
  inline friend DifNumber sec(const DifNumber &a)                     
  { return DifNumber(cos(a)).inverse(); }
  inline friend DifNumber asin(const DifNumber &a);
  inline friend DifNumber acos(const DifNumber &a);
  inline friend DifNumber atan(const DifNumber &a);

  inline friend DifNumber atan2(const DifNumber& y,const DifNumber& x) { return DifNumber(y).arcTangent(x); }
  inline friend DifNumber atan2(const DifNumber& y,const double& x) { return DifNumber(y).arcTangent(DifNumber(x)); }
  inline friend DifNumber atan2(const double& y,const DifNumber& x) { return DifNumber(y).arcTangent(x); }

  inline friend DifNumber exp(const DifNumber& a)                      
  { double e=exp(a.number()); return DifNumber(e,a,e); }
  inline friend DifNumber cosh(const DifNumber& a)                     
  { return 0.5*(DifNumber(exp(a))+=exp(-a)); }
  inline friend DifNumber sinh(const DifNumber& a)                     
  { return 0.5*(DifNumber(exp(a))-=exp(-a)); } 
  inline friend DifNumber tanh(const DifNumber& a)                     
  { double t=tanh(a.number()); return DifNumber(t,a,1.0-t*t); }
  // inline friend DifNumber sqrt(const DifNumber& a)                     
  //{ DifNumber temp(a); return temp.squareRoot(); }
  inline friend DifNumber sqrt(const DifNumber& a)                     
  { return DifNumber(a).squareRoot(); }
  inline friend DifNumber log(const DifNumber& a)                      
  { return DifNumber(log(a.number()),a,1.0/a.number()); }
  inline friend DifNumber fabs(const DifNumber& a)                     
  { return DifNumber(fabs(a.number()),a,a.sign()); }
  inline friend DifNumber pow(const DifNumber& a,const DifNumber& b)   { return DifNumber(a).power(b); }
  inline friend DifNumber pow(const DifNumber& a,int i)                { return pow(a,(double)i);}
  inline friend DifNumber pow(const DifNumber& a,float i)              { return pow(a,(double)i);}
  inline friend DifNumber pow(const DifNumber& a,double i)                {
                          return i==0?DifNumber(1.0,a,0.0)
                                   :DifNumber(pow(a.number(),i),a,i*pow(a.number(),i-1)); }


  // FIXME: this one should be inlined...
  friend double correlation(const DifNumber& a,const DifNumber& b,const HepSymMatrix& e);	// correlation for specified error
  
  inline friend double correlation(const DifNumber& a,const DifNumber& b); // correlation from default matrix
  

private:

  //data
  double _number;		// value of number
  int _npar;			// number of parameters
  const DifIndepPar* _indepPar;	// pointer to independent parameters 
  double _derivatives[_nmax];	// derivatives

  //functions
  inline DifNumber& copyDerivs(const DifNumber& n);  // copy derivatives from n
  inline DifNumber& scaleDerivs(const DifNumber& n,double factor); // copy byt scale by factor
  inline DifNumber& copyDerivs(const HepVector& v); // copy derivatives v
  inline DifNumber& check(const DifNumber& a); // check for consistent numbers
  inline DifNumber& setLike(const DifNumber& a); // set change to same type as |a|
};


//io
inline std::ostream& operator<<(std::ostream& o,const DifNumber& n) { n.print(o); return o; }

DifNumber solveQuad(const DifNumber& a, const DifNumber& b, const DifNumber& c, int pref, Code& code);

#include "difAlgebra/DifNumber.icc"

#endif
