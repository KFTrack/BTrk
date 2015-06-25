// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetSurfaceElem.cc,v 1.26 2002/12/30 15:44:29 dbrown Exp $
//
//  Description:
//  DetElem for 2-D surface objects
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//  
//  Authors: Dave Brown, 10/8/96
//------------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/ErrLog.hh"
#include "BTrk/DetectorModel/DetSurfaceElem.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/DetectorModel/Intersection.hh"
#include "BTrk/DetectorModel/DetMaterial.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <assert.h>
using namespace CLHEP;
//
static const double _epsilon = 1.0e-3;

DetSurfaceElem::DetSurfaceElem(DetSurfaceType* type,const char* oname,int id,
			       const DetSurface& surf) :
  _osurf(surf.copyOf()),DetElem(type,oname,id)
{
  _etrans = _osurf->transform();

}
//
//  Intersect uses the Intersection helper class
//
int
DetSurfaceElem::intersect(const Trajectory* trktraj,DetIntersection& dinter) const {
  int jflag = 0;
  SurfacePoint surfcoord;
  Intersection intersection( *trktraj, *_osurf );
  double flightrange[2];
  double flightdist = trktraj->lowRange();//
//  Use the initial values to define the starting point and the search range,
//  so long as that is valid
//
  if(dinter.pathrange[1] > dinter.pathrange[0]){
    flightdist = dinter.pathlen;
    flightrange[0] = dinter.pathrange[0];
    flightrange[1] = dinter.pathrange[1];
  } else {
    flightrange[0] = trktraj->lowRange();
    flightrange[1] = trktraj->hiRange();
  }
  TrkErrCode iflag = intersection.intersect( flightdist, surfcoord,trkOut,flightrange);
//
//  loop over successful intersections (the first one might not be the
//  right one!
//
  while(iflag.success()){
    TwoDCoord detcoord(surfcoord.array());
    jflag = _dtype->physicalMaterial(&detcoord);
    if(jflag){
//
//  Determine the flight range for traversing this element
//
//  This next section is temporary, and should be replaced by a more
//  robust interface for determining the entrance/exit of a track
//  going through a DetSurfaceType.
//
      double thickness = surfaceType()->thickness(&detcoord);
      Hep3Vector direction = trktraj->direction(flightdist);
      thickness /= fabs(direction.dot(_osurf->normal(surfcoord))); // correct for incident angle
      dinter = DetIntersection((DetElem*)this,trktraj,
				    flightdist,
				    flightdist-thickness/2,
				    flightdist+thickness/2);
      dinter.flag[0]=-2;
      dinter.flag[1]=-2;
//
      break;
    } else { // try again
      flightrange[0] = flightdist + _epsilon;
      iflag = intersection.intersect( flightdist, surfcoord,trkOut,flightrange);
    }
  }
  return jflag;
}
//
//  Outline functions
//

void
DetSurfaceElem::physicalOutline(std::vector<HepPoint>& pvec) const {
//
//  Reset the list
//
  pvec.clear();
//
//  Get the surface coordinate outline
//
  const std::vector< TypeCoord* >* outline = _dtype->outline();
//  Loop over the outline
  for(int icorn=0;icorn<outline->size();icorn++){
//  Turn them into space points, and add them to the vector
    pvec.push_back(_osurf->spacePoint((*outline)[icorn]->array()));
  }
}


//
const DetMaterial&
DetSurfaceElem::material(const DetIntersection& dinter) const {
//
//  Find the surface point
//
  HepPoint spacepoint = dinter.trajet->position(dinter.pathlen);
  SurfacePoint surfpoint;
  int iflag = _osurf->surfacePoint(spacepoint,surfpoint,0.2);
//
//  From the position on the wafer, ask the type how what material was traversed
//
  TwoDCoord detcoord(surfpoint.array());
  return detectorType()->material(&detcoord);
}

// Protected member functions
HepPoint
DetSurfaceElem::coordToPoint( const TypeCoord* aCoord ) const
{
  return _osurf->spacePoint( SurfacePoint(aCoord->array()));
}
