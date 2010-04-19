// $Id: BFieldFixed.hh 497 2010-01-14 09:06:53Z stroili $
// Description:	Class Header for |BFieldFixed|
//              Provide an arbitary fixed field
// Author List:A. Snyder, Copyright (C) 1998	SLAC
#ifndef BFieldFixed_HH
#define BFieldFixed_HH

#include <iostream>
#include "BaBar/BaBar.hh"
#include "BField/BField.hh"
#include "CLHEP/Vector/ThreeVector.h"

// class interface //
class BFieldFixed : public BField {

public:
  //construct
  BFieldFixed(double bx, double by, double bz,double delnom=0.0);
  BFieldFixed(const Hep3Vector& fieldvector,double delnom=0.0);

public:
  //destroy
  virtual ~BFieldFixed();

  //basics
  virtual Hep3Vector bFieldVect
  (const HepPoint &point=HepPoint(0,0,0))const;
  
  //print
  virtual void print(std::ostream &o)const;
// override the nominal field
  virtual double bFieldNominal()const;

private:
  Hep3Vector _field;
  double _delnom; // allow offsetting the nominal field
private:
  // Preempt copy constructor and operator=
  BFieldFixed&   operator= (const BFieldFixed&);
  BFieldFixed(const BFieldFixed &);
};
#endif





