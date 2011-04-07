// $Id: mu2eGradientField.hh 497 2010-01-14 09:06:53Z stroili $
// Description:	Class Header for |mu2eGradientField|
//              Provide an arbitary fixed field
// Author List:A. Snyder, Copyright (C) 1998	SLAC
#ifndef mu2eGradientField_HH
#define mu2eGradientField_HH

#include <iostream>
#include "BaBar/BaBar.hh"
#include "BField/BField.hh"
#include "CLHEP/Vector/ThreeVector.h"

// class interface //
class mu2eGradientField : public BField {

public:
// define the fields in the constant regions: the gradient is constant between those
  mu2eGradientField(double b0,double z0, double b1, double z1, double rmax);
public:
  //destroy
  virtual ~mu2eGradientField();

  //basics
  virtual Hep3Vector bFieldVect(const HepPoint &point=HepPoint(0,0,0))const;
protected:
private:
  Hep3Vector _b0,_b1;
  double _z0,_z1;
  double _rmax, _grad;
private:
  // Preempt copy constructor and operator=
  mu2eGradientField&   operator= (const mu2eGradientField&);
  mu2eGradientField(const mu2eGradientField &);
};
#endif





