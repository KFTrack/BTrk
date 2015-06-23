//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetElemSet.hh,v 1.4 1998/03/11 06:13:00 brownd Exp $
//
// Description:
//	DetElemSet Class - An Element Set is a DetSet which has 
//      a special relation with a single element, though it can
//      contain other elements and sets as any DetSet.  The
//      special element gives a spatial aspect to the set as a whole, and
//      permits an intermediate level of alignment between local and
//      global.  Intermediate alignment operates on all the contained
//      elements and sets (including the special element) as per global
//      alignment, taking the coordinate system of the special object as
//      the 'global' coordinate system.  These extra alignment functions
//      are the only difference between an DetElemSet and a generic Detector Set.
//      Note that this class can be used to give a spatial aspect to a
//      set without defining a prefered object simply by using a
//      DetNullElem object as the special object.
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
//      DNB       Jan 23, 1997  - convert SpatialSet to DetElemSet
//
// Copyright Information:
//	Copyright (C) 1997		Lawrence Berkeley Laboratory
//	Copyright (C) 1997	       CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------

#ifndef ELEMENTSET_HH
#define ELEMENTSET_HH

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Geometry/Transformation.h"
#include "DetectorModel/DetSet.hh"
class DetAlignElem;

//		---------------------
// 		-- Class Interface --
//		---------------------

class DetElemSet : public DetSet {

public:

// Constructors
  DetElemSet( const char*, int ,DetElem* );

// Destructor
  virtual ~DetElemSet();

// access
  const HepTransformation& transform() const { return specialelement->transform(); }
  const DetElem* specialElement() const { return specialelement; }
// 
  void applyIntermediate(const DetAlignElem&);
  void removeIntermediate(const DetAlignElem&);
//
private:
  DetElem* specialelement;
};

#endif

