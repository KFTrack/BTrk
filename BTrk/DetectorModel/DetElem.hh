// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetElem.hh,v 1.32 2004/12/14 07:10:17 bartoldu Exp $
//
//  Description:
//  Base Class to define a specific element of the tracking detector.  The types
//  must exist before elements can be created from them.  Also the subclasses
//  are responsable for providing the HepTransformation* to the base class
//  (note: the DetElem HepTransformation pointer MUST point to a valid
//  transform, which is owned (perhaps indirectly) by the element itself).
//  This class is the main interface between abstract
//  track representations (hits and fit parameters) and the 3-D material
//  model.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 7/17/96
// ------------------------------------------------------------------------------
#ifndef DETECTORELEMENT_HH
#define DETECTORELEMENT_HH
//----------------
// BaBar header --
//----------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/ErrLog.hh"
//
//  global includes
//
#include <iostream>
#include <string>
#include "BTrk/BbrGeom/HepPoint.h"
#include <vector>
//
//  Local includes
//
#include "BTrk/DetectorModel/DetType.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/TrkDirection.hh"
//
class TrkDifTraj;
class Trajectory;
class DetType;
class DetElem;
class DetAlignElem;
class HepTransformation;
class DetIntersection;
class GnuPlot;
//
//  Define the class
//
class DetElem{
public:
// 
//  Constructors
//
  DetElem();
  DetElem(const DetType*,const char*,int);
//  Assignment operator
  DetElem& operator = (const DetElem&);
//  equality operator (needed for old templates)
  bool operator == (const DetElem& other) const;
// see if an align elem matches this element
  bool match(const DetAlignElem& ae) const;
//
//  Destructor
//
  virtual ~DetElem();
//
//  Geometric functions
//
//  Intersect a trajectory with this detector element.  The DetIntersection
//  reference range limits initially define the search range and start value, but
//  are returned as the entrance and exit point of the intersection. If the input
//  lower limit is >= the upper limit, the range is ignored.
//
  virtual int intersect(const Trajectory*,DetIntersection&) const = 0;
// re-interesect with (potentially) an updated trajectory.  Default implementation just
// updates the trajectory in the intersection, but doesn't recompute any geometric quantities
  virtual bool reIntersect(const Trajectory*,DetIntersection&) const;
//
//  Material information from an intersection.  Default base class versions
//  are provided which assume any given intersection goes through only
//  one type of material.  If this is not the case, the element subclass
//  must overwrite these functions.  The user interface is only the second
//  function, which calls the first in the base class implementation.  This
//  allows subclasses which are not homogenous but still have only one
//  material type for a given intersection to use the base implementation
//  of materialInfo, only overwriting the base of material.
//
//  Note, this function now returns the (dimensionless) fractional
//  change in momentum associated with energy loss, not the energy
//  loss itself (ditto for the energy loss RMS).
//
// DNB 3/13/00 Added the particle direction as an argument to materialInfo;
// trkIn means the particle is passing inwards (pfrac>0),
// trkOut for passing outwards (pfrac<0) through the material.
public:
  virtual const DetMaterial& material(const DetIntersection&) const;
  virtual void materialInfo(const DetIntersection&,
			    double momentum,
			    TrkParticle const& tpart,
			    double& deflectRMS,
			    double& pFractionRMS,
			    double& pFraction,
			    trkDirection dedxdir=trkOut) const;
//
//  Alignment functions.  The name and ID number of the local DetAlignElem 
//  must match that of the DetElem.
//
  void applyGlobal(const DetAlignElem&);// apply global alignment
  void applyLocal(const DetAlignElem&); // apply local alignment
  void removeGlobal(const DetAlignElem&);// unapply global alignment
  void removeLocal(const DetAlignElem&); // unapply local alignment
  virtual void updateCache(); // update an elements cache (used for subclasses)
//
//  Access functions
//
  virtual void print(std::ostream& os) const;
  virtual void printAll(std::ostream& os ) const;
  const DetType* detectorType() const { return _dtype; }
  int elementNumber() const {return _ielem; }
  const std::string& elementName() const {return _ename; }
  const HepTransformation& transform() const { return *_etrans; }
//
//  Outline function
//
  virtual void physicalOutline(std::vector<HepPoint>&) const;
  virtual void gnuPlot( GnuPlot* ) const;
protected:
  virtual HepPoint coordToPoint( const TypeCoord* aCoord ) const = 0;
  // the ElemPointIterator class must be able to access coordToPoint function
  friend class DetElemPointIterator;  
  // Following is so derived classes can set transform through method
  HepTransformation*& myTransf() { return _etrans; }
  HepTransformation& transf() { return *_etrans; } // nonconst transform
  
  // private:
  int _ielem; // integer identifier; this is controled by the sub-class
  const std::string _ename; // name
  const DetType* _dtype; // pointer to the type
  HepTransformation* _etrans; // spatial transform
// this is used in the DetElemSet subclass of DetSet
  friend class DetElemSet;
};
#endif
