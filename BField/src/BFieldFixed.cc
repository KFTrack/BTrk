//--------------------------------------------------------------------------
//$Id: BFieldFixed.cc 497 2010-01-14 09:06:53Z stroili $
//Description: Class Implementation for |BFieldFixed|
//Author List:A. Snyder
//Copyright (C) 1998	SLAC
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "BField/BFieldFixed.hh"
using std::endl;
using std::ostream;


BFieldFixed::BFieldFixed(double bx, double by, double bz,double delnom):
  _field(bx,by,bz),_delnom(delnom)
{}

BFieldFixed::BFieldFixed(const Hep3Vector& fieldvector,double delnom):
  _field(fieldvector),_delnom(delnom)
{}

BFieldFixed::~BFieldFixed(){}

Hep3Vector
BFieldFixed::bFieldVect(const HepPoint &p)const
{
  return _field;
}

double
BFieldFixed::bFieldNominal()const {
  return _field.z() + _delnom;
}

void
BFieldFixed::print(ostream &o)const
{
  o << "BFieldFixed::" << endl;
  o << "Bx:" << _field.x() << endl;
  o << "By:" << _field.y() << endl;
  o << "Bz:" << _field.z() << endl;
  o << "Nominal field offset:" << _delnom << endl;
}

