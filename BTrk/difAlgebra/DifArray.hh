//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DifArray.hh 501 2010-01-14 12:46:50Z stroili $
//
// Description:
//	Class Header for |DifArray|
//      What i do ?
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

#ifndef DifArray_HH
#define DifArray_HH

#include <assert.h>
#include <stdlib.h>

#include <iosfwd>
class DifNumber; 
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"

class DifArray {

public:
  DifArray(int n,int npar=0);	// make an array of n elements
  DifArray			// copy from a vector
  (const CLHEP::HepVector& v,int npar=0);
  DifArray(const DifArray& a);	// make a copy
  ~DifArray();			// destroy array
  DifNumber& operator[](int i);	// subscripting (elem=i)
  DifNumber& operator()(int i);	// subscripting (elem=i-1)
  DifNumber fetch(int i)const;	// fetch elem=i
  int nElem()const {return _nElem;} // number of elements
  DifArray& operator=(const DifArray&); // assignment operator
  CLHEP::HepMatrix jacobian()const;	// return matrix of derivatives
  void zero(int npar=0);	// zero this array
  void print(std::ostream& o)const;	// print this array

private:

  //data
  int _nElem;			// number of elements
  DifNumber* _pointer;		// pointer to data

  //private functions
  void copy(const DifArray& a);	// copy a to *this
  void copy(const CLHEP::HepVector& a,int npar=0); //copy a clhep vector
  

};

//io
inline std::ostream& operator<<(std::ostream& o,const DifArray& a){
  a.print(o);
  return o;
}


#endif
