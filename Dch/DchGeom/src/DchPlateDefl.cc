//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchPlateDefl.cc 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchPlateDefl
//      
//      
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1998	INFN & Padova University
//
// Revision History:
//	20020405  M. Kelsey -- Bug fix -- define copy ctor.
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchPlateDefl.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <fstream>
#include <assert.h>
#include <string>

//-------------------------------
// Collaborating Calss Headers --
//-------------------------------
#include "ErrLogger/ErrLog.hh"
#include "DchGeomBase/DchWirePar.hh"
using std::fstream;
using std::ifstream;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DchPlateDefl::DchPlateDefl() :
  _deflCorr(40, DchWirePar())
{
}

DchPlateDefl::DchPlateDefl(const char* file) :
  _deflCorr()
{
  const char* tagname = "Deflection";
  bool tag = false;
  std::string fline;

  ifstream alfile(file);
  if (alfile.good()) {
    while (!tag) {
      alfile >> fline;
      tag = (fline.find(tagname) != std::string::npos);
    }
    assert(tag);
  }
  if (alfile.eof()) {
    ErrMsg(fatal) << " alignment file empty..." << endmsg;
  }

  int id;
  double dr1, dr2, dphi1, dphi2, dz1, dz2, sag, diameter;

  for (int i = 0; i < 40; i++) { // one correction per layer
    assert(!(alfile.eof()));
    // first alignement is the global alignement
    alfile >> id >> dr1 >> dphi1 >> dz1 >> dr2 >> dphi2 >> dz2 >> sag
        >> diameter;
    _deflCorr.push_back(DchWirePar(id, dr1, dphi1, dz1, dr2, dphi2, dz2, sag,
        diameter));
  }

}

DchPlateDefl::DchPlateDefl(const DchPlateDefl& other) :
  _deflCorr(other._deflCorr)
{
  ;
}

DchPlateDefl::DchPlateDefl(std::vector<DchWirePar>& corrVec) :
  _deflCorr(corrVec)
{
  ;
}

DchPlateDefl::~DchPlateDefl()
{
}

//-------------
// Operators --
//-------------
DchPlateDefl&
DchPlateDefl::operator=(const DchPlateDefl& other)
{
  if (&other == this) return *this;
  _deflCorr.clear();
  _deflCorr = other._deflCorr;

  return *this;
}

//-------------
// Selectors --
//-------------
const DchWirePar&
DchPlateDefl::wireCorr(int i) const
{
  // no assert necessary: [] operator always tests bounds
  return _deflCorr[i];
}

