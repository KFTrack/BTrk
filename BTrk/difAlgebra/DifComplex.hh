//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DifComplex.hh 501 2010-01-14 12:46:50Z stroili $
//
// Description:
//	Class Header for |DifComplex|
//      Represent a complex diffential number
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	A. Snyder
//
// Copyright Information:
//	Copyright (C) 1997	SLAC
//
//------------------------------------------------------------------------Dif

#ifndef DifComplex_HH
#define DifComplex_HH

#include <assert.h>
#include <stdlib.h>
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/difAlgebra/DifNumber.hh"

#include <iosfwd>



class DifComplex {

public:

  //constructors
  inline DifComplex();
  inline DifComplex(const DifNumber& r);
  inline DifComplex(double r);
  inline DifComplex(const DifNumber& r,const DifNumber& i);
  inline DifComplex(double r,double i);
  inline DifComplex(const DifComplex& c);

  //access 
  inline const DifNumber& real()const {return _r;}
  inline DifNumber& gimeReal(){return _r;}
  inline const DifNumber& imag()const {return _i;}
  inline DifNumber& gimeImag(){return _i;}

  //set
  void setReal(const DifNumber& r) {_r=r;}
  void setImag(const DifNumber& i) {_i=i;}

  //operators and such
  inline DifComplex& operator+=(const DifComplex& a);
  inline DifComplex& operator-=(const DifComplex& a);
  inline DifComplex& operator*=(const DifComplex& a);
  inline DifComplex& operator/=(const DifComplex& a);
  inline DifComplex& operator/=(const DifNumber& a);
  inline DifComplex& operator/=(const double& a);
  inline DifComplex& operator=(const DifComplex& a);
  inline DifComplex operator-()const;

  //manupulate *this guy - reverse polish style
  inline DifComplex& conjugate();		// take complex conjugate
  inline DifComplex& inverse();		// invert
  inline DifComplex& flipsign();		// flip sign
  inline DifComplex& takeLog();		// take log
  inline DifComplex& power(const DifNumber& p); // raise to a power
  inline DifComplex& power(const DifComplex& p); // raise to a power

  //complex specfic function
  inline DifNumber mag()const;		// magnitude
  inline DifNumber magsq()const;	// magnitude squared
  inline DifNumber arg()const;		// phase
  inline DifNumber phase()const {return arg();} // a rose by ...

  //io
  inline void print(std::ostream& o)const;

  //friends
  inline friend DifComplex operator+(const DifComplex& a,const DifComplex& b);
  inline friend DifComplex operator-(const DifComplex& a,const DifComplex& b);
  inline friend DifComplex operator*(const DifComplex& a,const DifComplex& b);
  inline friend DifComplex operator/(const DifComplex& a,const DifComplex& b);
  inline friend bool operator==(const DifComplex& a,const DifComplex& b);
  inline friend bool operator!=(const DifComplex& a,const DifComplex& b);
  inline friend DifComplex sin(const DifComplex& a);
  inline friend DifComplex cos(const DifComplex& a);
  inline friend DifComplex tan(const DifComplex &a);
  inline friend DifComplex sec(const DifComplex &a);
  inline friend DifComplex asin(const DifComplex &a);
  inline friend DifComplex acos(const DifComplex &a);
  inline friend DifComplex atan(const DifComplex &a);
  inline friend DifComplex atan2(const DifComplex& a,const DifComplex& b);
  inline friend DifComplex exp(const DifComplex& a);
  inline friend DifComplex cosh(const DifComplex& a);
  inline friend DifComplex sinh(const DifComplex& a);
  inline friend DifComplex tanh(const DifComplex& a);
  inline friend DifComplex sqrt(const DifComplex& a);
  inline friend DifComplex log(const DifComplex& a);
  inline friend DifComplex pow(const DifComplex& a,const DifComplex& p);
  inline friend DifComplex pow(const DifComplex& a,const DifNumber& p);
  inline friend DifComplex cc(const DifComplex& a);
  inline friend DifComplex fromPolar(const DifNumber& mag,const DifNumber& phase);



  //destructor
  inline virtual ~DifComplex(){};

private:

  //data
  DifNumber _r;
  DifNumber _i;

  //functions
  inline DifNumber& gimeR() {return _r;}
  inline DifNumber& gimeI() {return _i;}
  inline const DifNumber& R()const {return _r;}
  inline const DifNumber& I()const {return _i;}
};

#include "BTrk/difAlgebra/DifComplex.icc"

inline std::ostream& operator<<(std::ostream& o,const DifComplex& c){
  c.print(o);
  return o;
}
 
#endif





