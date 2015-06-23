//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetElemSet.cc,v 1.9 2004/12/14 07:10:18 bartoldu Exp $
//
// Description:
//	Class DetElemSet
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//      David Brown, LBL
//
// History (add to end):
//      Gautier   Jan 14, 1997  - creation
//      DNB       Jan 23, 1997  - convert SpatialSet
//
// Copyright Information:
//	Copyright (C) 1997		Lawrence Berkeley Laboratory
//	Copyright (C) 1997	       CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include <assert.h>

#include "ErrLogger/ErrLog.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetElemSet.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Geometry/AngleSets.h"
#include "DetectorModel/DetAlignElem.hh"

//----------------
// Constructor  --
//----------------
DetElemSet::DetElemSet(const char* name, int sect, 
			DetElem* elem) 
  :   DetSet(name,sect),specialelement(elem)
{
  assert ( specialelement!=0 ); // require the element to exist
  *this += *specialelement;  // add the special element to the element list
}
//--------------
// Destructor --
//--------------
DetElemSet::~DetElemSet()
{;}

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

void
DetElemSet::applyIntermediate(const DetAlignElem& align) {
// This method applies the transformation
  HepTransformation transf(transform());
  HepTransformation Transwrtbabar(transform());
  transf *= align.transform();
  transf *= Transwrtbabar.inverse();
  //  transf *= transf.inverse();
  DetElemList allElemList;
  listAllElements(allElemList);
  DetElemList::const_iterator eiter = allElemList.begin();
  DetElem* elem = 0;
  while( eiter != allElemList.end() ) {
    elem = *eiter++;
    elem->transf().transform(transf);
    elem->updateCache();
  }
  _ready = false;  // set not ready : set elements have been modified
}

void
DetElemSet::removeIntermediate(const DetAlignElem& align) {
  HepTransformation transf(transform());
  HepTransformation Transwrtbabar(transform());
  transf *= align.inverseTransform();
  transf *= Transwrtbabar.inverse();
  //   transf *= transf.inverse();
  DetElemList allElemList;
  listAllElements(allElemList);
  DetElemList::const_iterator eiter = allElemList.begin();
  DetElem* elem = 0;
  while( eiter != allElemList.end() ) {
    elem = *eiter++;
    elem->transf().transform(transf);
    elem->updateCache();
  }
  _ready = false; // set not ready : set elements have been modified
}
