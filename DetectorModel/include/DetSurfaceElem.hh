// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetSurfaceElem.hh,v 1.16 2002/12/30 15:44:29 dbrown Exp $
//
//  Description:
//  DetElem for generic 2-D objects described by a surface.  The DetType pointed
//  to by the planar object must be subclassed to DetSurfaceType
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//  
//  Authors: Dave Brown, 10/8/96
//------------------------------------------------------------------------------
#ifndef SURFACEELEMENT_HH
#define SURFACEELEMENT_HH
//
//  Includes
//
#include <vector>
#include "CLHEP/Geometry/HepPoint.h"
#include "DetectorModel/DetElem.hh"
#include "DetectorModel/DetSurface.hh"
#include "DetectorModel/DetSurfaceType.hh"
//
//  Define the class
//
class DetSurfaceElem : public DetElem {
public:
//
//  Constructor
//
  DetSurfaceElem(DetSurfaceType*,const char*,int,const DetSurface&);
//  Destructor
  virtual ~DetSurfaceElem(){delete _osurf;}
//
//  Real versions of the virtual functions, see the base class for a description
//
  int intersect(const Trajectory*,DetIntersection&) const;
//
//  Access
//
  void physicalOutline(std::vector<HepPoint>&) const;
  virtual const DetSurface* surface() const { return _osurf; }
  HepPoint spacePoint( const TypeCoord& aCoord ) const {
    return coordToPoint(&aCoord); }
//
//  This next allows for non-homogenous surface elements
//
  const DetMaterial& material(const DetIntersection&) const;
protected:
  virtual HepPoint coordToPoint( const TypeCoord* aCoord ) const;
  DetSurfaceType* surfaceType() const {
    return (DetSurfaceType*)detectorType(); }
//
private:
  DetSurface* _osurf; // object surface.
};
#endif
