// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetNullElem.cc,v 1.14 2002/12/30 15:44:28 dbrown Exp $
//
//  Description:
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//  
//  Authors: Gautier Hamel de Monchenault, 1/15/97
//------------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/ErrLog.hh"
#include <assert.h>
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/DetectorModel/DetNullElem.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/BbrGeom/Transformation.h"
//

DetNullElem::DetNullElem(const char* name,int id,
			       const HepTransformation& transf) :
  DetElem(0,name,id)
{
  _etrans = new HepTransformation(transf);
}

DetNullElem::~DetNullElem() 
{
  delete _etrans;
}
  
HepPoint 
DetNullElem::coordToPoint( const TypeCoord* ) const 
{ 
  return HepPoint(); 
}



//  Intersect uses the Intersection helper class - trivial implementation
int
DetNullElem::intersect(const Trajectory* trktraj,DetIntersection& dinter) const {
  return 0;
}

//  Outline
void
DetNullElem::physicalOutline(std::vector<HepPoint>& stlvec) const {
  stlvec.clear();
}

