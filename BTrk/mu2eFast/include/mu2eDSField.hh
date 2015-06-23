// $Id: mu2eDSField.hh 497 2010-01-14 09:06:53Z stroili $
// Description:	Class Header for |mu2eDSField|
//              Provide an arbitary fixed field
// Author List:A. Snyder, Copyright (C) 1998	SLAC
#ifndef mu2eDSField_HH
#define mu2eDSField_HH

#include <iostream>
#include "BaBar/BaBar.hh"
#include "BField/BField.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "mu2eFast/BFMap.hh"

// class interface //
class mu2eDSField : public BField {

public:
  //construct from file; optionally amplify the distortions by the given factor
  mu2eDSField(const std::string& mapfile,double dfactor=1.0);
public:
  //destroy
  virtual ~mu2eDSField();

  //basics
  virtual Hep3Vector bFieldVect
  (const HepPoint &point=HepPoint(0,0,0))const;
  // override the nominal field
  virtual double bFieldNominal()const;
protected:
  // Read a MECO GMC format map.
  void readGMCMap( const std::string& filename, mu2e::BFMap& bfmap );

  // Compute the size of the array needed to hold the raw data of the field map.
  
  int computeArraySize( int fd, const std::string& filename );

private:
  
  // coordinate conversions
  static Hep3Vector trackerCenterInMu2eCoordinates;
  static Hep3Vector trackerCenterInFastSimCoordinates;
  // the actual field
  mu2e::BFMap _fieldmap;
  double _bnom;
  double _dfactor;
  bool _distort;
private:
  // Preempt copy constructor and operator=
  mu2eDSField&   operator= (const mu2eDSField&);
  mu2eDSField(const mu2eDSField &);
};
#endif





