//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchTestGeom.hh 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//      does some tests on Dch geometry
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: R. Stroili
//
// Copyright (C)  1999  INFN 
// 
//------------------------------------------------------------------------

#ifndef DCHTESTGEOM_HH
#define DCHTESTGEOM_HH

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"
#include "Framework/AbsParmBool.hh"
#include "Framework/AbsParmDouble.hh"
#include "Framework/AbsParmGeneral.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"

class DchDetector;
class HepTupleManager;
class HepHistogram;
class HepTuple;

// Class interface //
class DchTestGeom : public AppModule {

public:
  DchTestGeom(const char* const theName, const char* const theDescription );
  virtual ~DchTestGeom();

  AppResult           beginJob( AbsEvent* anEvent );
  AppResult           event   ( AbsEvent* anEvent );
  AppResult           endJob  ( AbsEvent* anEvent );

  void                RemedyWoe(double initMom, double phi);

private:	
  const DchDetector* _dch;
  
protected:
  // local parameter variables
  AbsParmBool      _dump;
  AbsParmBool      _printwires;
  AbsParmBool _intersectdch; // run intersection test of SvtDetector
  AbsParmGeneral<int> _nct; // number of cos theta scan points
  AbsParmGeneral<int> _np; // number of phi scan points
  AbsParmBool _validatenext; // validate nextIntersection
  AbsParmDouble _ctlow; // cos theta scan lower limit
  AbsParmDouble _cthi; // cos theta scan upperlimit
  AbsParmDouble _plow; // phi scan lower limit
  AbsParmDouble _phi; // phi scan upper limit
  AbsParmDouble _tlen; // trajectory length
  AbsParmGeneral<HepString> _trajtype; // type of trajectory (either line or helix)
  AbsParmGeneral<HepString> _set; // set to test (gasset or layerset)
  AbsParmDouble _curve; // curvature; vaule 0.0 means line trajectory
  AbsParmDouble _nphi; // # of phi slices for David's scan 
  AbsParmDouble _mom; // momentum for David's scan 
  AbsParmBool _scanVolume; // David's scan volume function
  HepTuple* _intersectntup; // intersection test results ntuple

//   const DchDetector* _dch;
};
#endif







