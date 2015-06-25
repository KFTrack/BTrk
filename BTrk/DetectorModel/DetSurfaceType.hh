// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetSurfaceType.hh,v 1.7 2002/12/30 15:44:29 dbrown Exp $
//
//  Description:
//  DetType subclass for surface objects
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 1/7/97
//------------------------------------------------------------------------------
#ifndef DETSURFACETYPE_HH
#define DETSURFACETYPE_HH
#include "BTrk/DetectorModel/DetType.hh"
#include <vector>

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class TwoDCoord;

//
//  Define the class
//
class DetSurfaceType : public DetType {
public:
  DetSurfaceType();
  DetSurfaceType(const char*,int);
  //
  //  New constructor with a list of corners
  //
  DetSurfaceType(const char*, int, std::vector< TwoDCoord >& );

  virtual ~DetSurfaceType();
//
//  Thickness of the type normal to the surface coordinate
//
  virtual double thickness(const TwoDCoord* point) const = 0;
//
//  Effective thickness given the direction cosines describing the
//  incident angle along the two (orthogonal) coordinate directions 
//
  virtual double effectiveThickness(const TwoDCoord* point,
				    double dircos1,
				    double dircos2) const;

protected:

  //
  //  Returns true if thisPoint is inside physical limits.
  //  Warning : the function is not pure virtual to prevent
  //  users for having to implement it in all derivatives.
  //  The default implementation is not general, though.
  //  It assumes that the boundary is made of segments, 
  //  determined by the points stored in pointStore, 
  //  that the origin is inside physical limits and that
  //  angles between two consecutive segments are of constant sign.
  //  If your surface type doesn't obey these criteria, 
  //  the function must be overwritten.
  //  A general base implementation will be provided eventually.
  //
  virtual bool insideLimitsOf( const TwoDCoord& thisPoint ) const;

  // Helper functions
  virtual bool insideLine( const TwoDCoord& thisPoint,
			   const TwoDCoord& p1,
			   const TwoDCoord& p2 ) const;
  
private:

};
#endif
