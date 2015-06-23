//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: GenDchGeom.hh 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//     Creates a DchGeom object and stores it in AbsEnv
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//
// Copyright (C)  1996  The Board of Trustees of  
// The Leland Stanford Junior University.  All Rights Reserved.
//------------------------------------------------------------------------

#ifndef GENDCHGEOM_HH
#define GENDCHGEOM_HH

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"
#include "Framework/AbsParmFilename.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"

class DchDetector;
class HepTupleManager;
class HepTuple;

// Class interface //
class GenDchGeom : public AppModule {

public:
  GenDchGeom(const char* const theName, const char* const theDescription );
  virtual ~GenDchGeom();

  AppResult           beginJob( AbsEvent* anEvent );
  AppResult           event   ( AbsEvent* anEvent );
  AppResult           endJob  ( AbsEvent* anEvent );

private:	
  HepTupleManager* _hFile;
  HepTuple*        _ntple;
  
protected:
  // local parameter variables
  AbsParmFilename _gfile;
  AbsParmIfdStrKey _trkRecoTrkList;
//   const DchDetector* _dch;
};
#endif







