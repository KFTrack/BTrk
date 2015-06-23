// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetNullElem.hh,v 1.13 2002/12/30 15:44:28 dbrown Exp $
//
//  Description:  DetNullElem is a minimal concrete subclass
//                of DetElem, with trivial implementations
//                for the tracking methods, and simple implementations
//                for the outline methods. 
//                DetNullElem owns a HepTransformation, the minimum
//                requirement to be a DetElem.
//
// Author List:
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//
// History (add to end):
//      Gautier   Jan 14, 1997  - creation
//      DNB       Jan 24, 1997 - make DetNullElem own it's transform, change constructor
//
// Copyright Information:
//	Copyright (C) 1997		Lawrence Berkeley Laboratory
//	Copyright (C) 1997	       CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------
#ifndef NULLELEMENT_HH
#define NULLELEMENT_HH

#include "CLHEP/Geometry/HepPoint.h"
#include "DetectorModel/DetElem.hh"
//
//  Class definition
//
class DetNullElem : public DetElem {
public:

//  Constructor

  DetNullElem(const char*, int, const HepTransformation&);

//  Destructor
  virtual ~DetNullElem();

//  Real versions of the virtual functions
  int intersect(const Trajectory*,DetIntersection&) const;
//  Access
  virtual void physicalOutline(std::vector<HepPoint>&) const;
private:
  HepPoint coordToPoint( const TypeCoord* ) const;
// prohibit
  DetNullElem& operator = (const DetNullElem& other);
  DetNullElem(const DetNullElem&);
};
#endif
