//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DifRotation.hh 501 2010-01-14 12:46:50Z stroili $
//
// Description:
//	Class Header for |DifRotation|
//      Rotation matrix and rotate things
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

#ifndef DifRotation_HH
#define DifRotation_HH

#include "BTrk/difAlgebra/DifVector.hh"

class DifRotation {

public:

  //constructors
  DifRotation();		// unit matrix 
  DifRotation			// construct from euler
  (const DifNumber& alpha,
   const DifNumber& beta,
   const DifNumber& gamma);
  DifRotation(double, double, double);

  DifRotation			// constructionr from axes
  (const DifVector& xp,const DifVector& yp,const DifVector& zp);
  DifRotation			// z=x X y
  (const DifVector& xp,const DifVector& yp);

  ~DifRotation() {};

  //access
  inline DifVector xnew()const {return _xnew;}
  inline DifVector ynew()const {return _ynew;}
  inline DifVector znew()const {return _znew;}

  //rotate a vector
  void rotate(DifVector& v)const;

  //error check
  int fail()const;		// check for orthonormality failure

private:

  //data - store as vector
  DifVector _xnew;
  DifVector _ynew;
  DifVector _znew;
};

#endif
