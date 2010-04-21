// ---------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetSurfaceSet.cc,v 1.53 2004/12/14 07:10:18 bartoldu Exp $
//
//  Description:
//  Special form of a DetSet where the elements are arranged roughly around
//  a surface.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/1/97
//---------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "ErrLogger/ErrLog.hh"
#include "BaBar/Constants.hh"
#include <math.h>
#include <cfloat>
#include <assert.h>
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetSurfaceSet.hh"
#include "DetectorModel/DetAlignElem.hh"
#include "DetectorModel/DetElem.hh"
#include "BbrGeom/Trajectory.hh"
#include "DetectorModel/Intersection.hh"
#include "CLHEP/Geometry/Transformation.h"

#include <vector>
#include <algorithm>
using std::endl;
using std::ostream;

static const double _epsilon = 1.0e-6;
static const double _surfdelta = 0.01; // roundoff protection
static const double _perpdelta = 0.01; // roundoff protection
//
//  Constructor just sets up the storage and lists.  The hash table organizes the
//  contained elements into an evenly spaced 2-d grid.
//
DetSurfaceSet::DetSurfaceSet(const char* name,int IDNum,
			     const DetSurface& seed,
			     int nxHash,int nyHash) :
  DetSet(name,IDNum),_refsurface(seed.copyOf()),_nhash(nxHash*nyHash),
//PE  _hash(0)
  _theHash(nxHash*nyHash)
{
//
//  Set single coordinate info
//
  _nbuckets[0] = nxHash;
  _nbuckets[1] = nyHash;
  _distlimits[0] = -1.0;
  _distlimits[1] = 1.0;
  _surfaces[0] = 0;
  _surfaces[1] = 0;
  for(int icoord=0;icoord<2;icoord++){
    _wrapped[icoord] = _refsurface->wrappedCoordinate(icoord);
    _nwrap[icoord] = 0;
    _surfaces[icoord] = 0;
    _limits[0][icoord] = FLT_MAX;
    _limits[1][icoord] = -FLT_MAX;
    _range[icoord] = 0.0;
  }
}
// Destructor; clean up the surfaces and hash table
DetSurfaceSet::~DetSurfaceSet(){
  delete _surfaces[0];
  delete _surfaces[1];
  delete _refsurface;
}
//
//  implementations for testing readiness 
//
bool
DetSurfaceSet::isReady() const {
  if(!_ready) return false;
// test all the subsets recursively
  if(_slist.size() > 0) {
    DetSetList::const_iterator siter = slist().begin();
    DetSet* set = 0;
    while( siter != slist().end() ) {
      set = *siter++;
      if(!set->isReady()) return false;
    }
  }
  return true;
}

void
DetSurfaceSet::makeReady() const {
  if(!_ready){
    buildHash();    // create the hash table
    ready() = true;
//  make the subsets ready recursively
    if(_slist.size() > 0) {
      DetSetList::const_iterator siter = slist().begin();
      DetSet* set = 0;
      while( siter != slist().end() ) {
	set = *siter++;
	set->makeReady();
      }
    }
  }
}
//
//  Intersection function.  This first determines an approximate intersection point
//  using the surfaces.  The trajectory is then mapped onto the surface over its
//  'thickness', making a 2-D line segement (this is to insure that we don't loose
//  any intersected elements due to the non-cylindrical nature of individual elements).
//  The segment is then hashed to find the set of likely intersected elements (note
//  that this pre-selection is very conservative, and should never loose any elements).
//  The sub-set of likely elements is tested for uniqueness, then the explicit
//  intersection is called for each.
//

// The same with an STL container in the interface

void
DetSurfaceSet::intersection(std::vector<DetIntersection>& divec,
			    const Trajectory* traj,double* myrange,
			    bool clear) const {
//
//  clear the list at the top level call
//
  if(clear)divec.clear();
//
//
//  First, make sure the hash table is up-to-date
//
  if(!isReady()) makeReady();
//

//  Setup the range 
//
  double range[2];
    if(myrange != 0){
    range[0] = myrange[0];
    range[1] = myrange[1];
  } else {
    range[0] = traj->lowRange();
    range[1] = traj->hiRange();
  }
//
//  Intersect the trajectory with the surfaces
//
  double flightdist[2];
  flightdist[0] = flightdist[1] = range[0];
  SurfacePoint surfinter[2];
  double distrange[2];
  bool between(false);
  bool hitsurfs = surfaceIntersections(traj,surfinter,flightdist,distrange,
				       range,between);
//
//  Loop while the surfaces give intersections, to cover the case of loopers
//
  while(hitsurfs){
//
//  Require at least 1 point within the surface range of interest
//
    if(withinRange(surfinter[0]) ||
       withinRange(surfinter[1]) ){
//
//  Set the range to be between the surfaces (up to the trajectory limits)
//
      double finalrange[2];
      finalrange[0] = std::min(flightdist[0],flightdist[1]);
      finalrange[1] = std::min(std::max(flightdist[0],flightdist[1]),range[1]);
//
//  Build the list of element candidates from the cells
//
      // intelements needs to be a SORTED vector
      static std::vector<const DetElem*> intelements(16);
      interElements(surfinter,intelements,distrange);

//
//  Actually intersect the elements on the reduced list, and insert it if needed
//
      for ( std::vector<const DetElem*>::const_iterator elemIter = intelements.begin();
	    elemIter != intelements.end();
	    ++elemIter ) {
	const DetElem* elem = *elemIter;
	DetIntersection dinter(0,traj,finalrange);
//
//  Loop over intersections, for the (rare) case of glancing intersections
//  on a single element
//
	int nelem = 0;
	while(elem->intersect(traj,dinter)){
//
//  Check that the element intersection lies in the correct range
//
	  if(dinter.pathlen>= finalrange[0] &&
	     dinter.pathlen<= finalrange[1]){
	    divec.push_back(dinter);
	    nelem++;
	    if(nelem > 10){
	      ErrMsg(error) << "DetSurfaceSet: element iteration truncated" << endmsg;
	      break;
	    }
	    dinter.pathrange[0] = dinter.pathrange[1] + _epsilon;
	    dinter.pathrange[1] = std::max(dinter.pathrange[0]+_epsilon,finalrange[1]);
	  } else
	    break;
//
//  Check if we're outside the range
//
	  if(dinter.pathrange[0] > finalrange[1])break;
	}
      }
    }

    std::sort(divec.begin(), divec.end());

//
//  Try re-intersecting the surfaces, to catch the loopers
//
    double frange[2];
    frange[0] = std::max(flightdist[0],flightdist[1]) + _epsilon;
    flightdist[0] = flightdist[1] = frange[0];
    frange[1] = range[1];
    hitsurfs = (frange[0]<frange[1]) && surfaceIntersections(traj,surfinter,
							     flightdist,distrange,
							     frange,between);
  }
}


//
//  Find the first intersection within the specified range.
//
bool
DetSurfaceSet::firstIntersection(const Trajectory* traj,
				 DetIntersection& first,
				 double* myrange) const {
//
//  Setup the range
//
  double frange[2];
  if(myrange != 0){
    frange[0] = myrange[0];
    frange[1] = myrange[1];
  } else {
    frange[0] = traj->lowRange();
    frange[1] = traj->hiRange();
  }
//
//  Intersect the trajectory with the surfaces.
//
  double flightdist[2];
  bool between;
  flightdist[0] = flightdist[1] = frange[0];
  SurfacePoint surfinter[2];
  double distrange[2];
  bool hitsurfs = (frange[0]<frange[1]) && 
    surfaceIntersections(traj,surfinter,flightdist,distrange,frange,between);
  bool foundelem = false;
  if(hitsurfs){
//
//  Require at least 1 point within the surface range of interest
//
    if(withinRange(surfinter[0]) ||
       withinRange(surfinter[1]) ){
//
//  Set the range to be between the surfaces (up to the trajectory limits)
//
      double minpath = FLT_MAX;
      double finalrange[2];
      finalrange[0] = std::max(std::min(flightdist[0],flightdist[1]),frange[0]);
      finalrange[1] = std::min(std::max(flightdist[0],flightdist[1]),frange[1]);
//
//  Build the list of element candidates from the cells
//
      // intelements needs to be a SORTED vector.
      static std::vector<const DetElem*> intelements(16);
      interElements(surfinter,intelements,distrange);
//
//  Actually intersect the elements on the reduced list
//
      for ( std::vector<const DetElem*>::const_iterator elemIter = intelements.begin();
	    elemIter != intelements.end();
	    ++elemIter ) {
	const DetElem* elem = *elemIter;
	DetIntersection dinter(0,traj,finalrange);
	if(elem->intersect(traj,dinter))
	  if( dinter.pathlen < minpath){
	    minpath = dinter.pathlen;
// further restrict the range: we now only care about intersections
// _closer_ than the one we just found
	    finalrange[1] = minpath;
	    first = dinter;
	    foundelem = true;
	  }
      }
    }
  }
  return foundelem;
}
//
//  Setup the surfaces and hash tables
//
void
DetSurfaceSet::buildHash() const {
  
  _theHash.clear();
  for(int ihash=0;ihash<_nhash;ihash++) {
    _theHash.push_back( std::vector<DetElemSurfRange>());
  }

  if(_surfaces[0])delete surfaces(0);
  if(_surfaces[1])delete surfaces(1);
//
//  Locate all the elements in this and all contained sets
//
  DetElemList allElemList;
  listAllElements(allElemList);
//
//  Find the global limits
//
  if(allElemList.size() > 0)
    globalLimits(allElemList);
//
//  Build the surfaces which enclose the elements.
//
  for(int isurf=0;isurf<2;isurf++) {
    surfaces(isurf) = _refsurface->copyOf(_distlimits[isurf]); }
// reset the reference to be centered
  double midpoint = (_distlimits[0] + _distlimits[1])/2.0;
  DetSurface* newref = _refsurface->copyOf(midpoint);
  delete refSurface();
  refSurface() = newref;
  distlimits(0) -= midpoint;
  distlimits(1) -= midpoint;
//
//  recalculate the global limits again after modification of reference surface
//
  if(allElemList.size() > 0)
    globalLimits(allElemList);
//
//  Loop over the elements and insert them in the hash table
//
  DetElemList::const_iterator iter = allElemList.begin();
  DetElem* elem = 0;
  while( iter != allElemList.end() ) {
    elem = *iter++;
//
//  Add this element to the hash table
//
    hashAddElem(elem);
  }

  // Since switching to STL, do a sort on all the contained vectors,
  // to agree with previous convention.
  std::vector<std::vector<DetElemSurfRange> >::iterator hiter2 = _theHash.begin();
  while( hiter2 != _theHash.end() ) {
    std::vector<DetElemSurfRange>& sortVec = *hiter2++;
    std::sort( sortVec.begin(), sortVec.end() );
  }
}
//
//  check the range
//
bool
DetSurfaceSet::withinRange(const SurfacePoint& spoint) const {
  bool within = true;
  double coordvalue;
  for(int icoord=0;icoord<2;icoord++){
    coordvalue = spoint[icoord];
    if(_wrapped[icoord])
//
//  Try to push the point into range
//
      if(coordvalue < _limits[0][icoord])
	coordvalue += Constants::twoPi;
      else if(coordvalue > _limits[1][icoord])
	coordvalue -= Constants::twoPi;
    within &= coordvalue >= _limits[0][icoord] &&
      coordvalue <= _limits[1][icoord];
  }
  return within;
}
//
//  alignment functions
//
void
DetSurfaceSet::applyGlobal(const DetAlignElem& glob){
//
//  Fall back to the default version for the elements
//
  if(!isReady()) makeReady();
  DetSet::applyGlobal(glob);
//
//  Update the surfaces if they exist
//
  for(int isurf=0;isurf<2;isurf++)
    if(_surfaces[isurf])
      _surfaces[isurf]->transform()->transform(glob.transform());
}
//
void
DetSurfaceSet::removeGlobal(const DetAlignElem& glob){
  if(!isReady()) makeReady();
  DetSet::removeGlobal(glob);
  for(int isurf=0;isurf<2;isurf++)
    if(_surfaces[isurf])
      _surfaces[isurf]->transform()->transform(glob.inverseTransform());
}
//
//  Given a (spatial) perimeter, translate it into a SurfaceCoordinate
//  range limit (and average.  This takes into account wrapped coordinates
//  by making sure that the outline is continuous in surface space.  This
//  also returns the maximum and minimum perpendicular distance from the
//  surface to the perimeter.
//
void
DetSurfaceSet::elemSurfaceOutline(std::vector<HepPoint>& hpvec,
				  SurfacePoint& minsurf,
				  SurfacePoint& maxsurf,
				  SurfacePoint& midsurf,
				  double& mindist,
				  double& maxdist ) const {
//
//  Initialize the limits
//
  mindist = FLT_MAX;
  maxdist = -FLT_MAX;
  int icoord;
  for(icoord=0;icoord<2;icoord++){
    minsurf[icoord] = FLT_MAX;
    maxsurf[icoord] = -FLT_MAX;
  }
//
//  Convert the first point to surface coordinates; this initializes
//  the 'reference point' for getting the coordinate wrapping right.
//
  double sdist;
  SurfacePoint refsurf,surf;
  Hep3Vector norm;
  sdist = project(hpvec[0],norm,refsurf);
//
//  Loop over the points in pairs (edge segments)
//
  int npoint = hpvec.size();
  for(int ipoint=0;ipoint<npoint;ipoint++){
    int jpoint = (ipoint+1) % npoint; // next point
//
//  Get the perp limits for this segment, and update the global limits
//
    double segmin,segmax;
    _refsurface->segmentMinMax(hpvec[ipoint],hpvec[jpoint],
			       segmin,segmax);
    mindist = std::min(mindist,segmin);
    maxdist = std::max(maxdist,segmax);
//
//  Get the surface coordinate of this point
//
    sdist = project(hpvec[ipoint],norm,surf);
//
//  Get the wrapping right.  This assumes that adjacent points
//  will always be on the same side of 2pi.
//
    for(icoord=0;icoord<2;icoord++){
      if(_wrapped[icoord]){
	if(fabs(surf[icoord]-refsurf[icoord]) > Constants::pi){
	  if(surf[icoord]-refsurf[icoord] > Constants::pi)
	    surf[icoord] -= Constants::twoPi;
	  else
	    surf[icoord] += Constants::twoPi;
	}
      }
//
//  Update the limts
//
      minsurf[icoord] = std::min(minsurf[icoord],surf[icoord]);
      maxsurf[icoord] = std::max(maxsurf[icoord],surf[icoord]);
    }
//
//  Set the reference for the next point
//
    refsurf = surf;
  }
//
//  Find the midpoint
//
  for(icoord=0;icoord<2;icoord++)
    midsurf[icoord] =( minsurf[icoord]+maxsurf[icoord])/2.0;
}
//
//  Find the global limits
//
void
DetSurfaceSet::globalLimits( DetElemList& allElemList) const {
//
//  Initialize limits
//
  int icoord;
  for(icoord=0;icoord<2;icoord++){
    limits(0,icoord) = FLT_MAX;
    limits(1,icoord) = -FLT_MAX;
  }
  distlimits(0) = FLT_MAX;
  distlimits(1) = -FLT_MAX;
//
//  Loop over the elements (including the elements all of the subsets),
//  and find the global limits (and in and out of the surface)
//
  DetElemList::const_iterator iter = allElemList.begin();
  std::vector<HepPoint> hpvec;
  double mindist,maxdist;
  SurfacePoint minsurf,maxsurf,midsurf,refsurf;
//
//  Get the surface point for the first element, as a reference for
//  resolving wrapping ambiguitites.  Protect against the first element
//  being a null element.
//
  DetElem* elem = *iter++;
  elem->physicalOutline(hpvec);
  while(hpvec.size() == 0){
    elem = *iter++;
    elem->physicalOutline(hpvec);
  }
  elemSurfaceOutline(hpvec,minsurf,maxsurf,refsurf,mindist,maxdist);
//
//  Loop over the elements
//
  iter = allElemList.begin();
  while( iter != allElemList.end() ) {
    elem = *iter++;
//
//  Get the physicalOutline perimeter of the element
//
    elem->physicalOutline(hpvec);
//
//  Convert this to limits for this element (if it's not a null element)
//
    if(hpvec.size() > 0){
      elemSurfaceOutline(hpvec,minsurf,maxsurf,midsurf,mindist,maxdist);
//
//  Resolve any wrapping ambiguities for wrapped coordinates
//
      for(icoord=0;icoord<2;icoord++){
	if(_wrapped[icoord]){
	  if(fabs(midsurf[icoord]-refsurf[icoord])>Constants::pi){
	    if(midsurf[icoord]-refsurf[icoord] > Constants::pi){
	      minsurf[icoord] -= Constants::twoPi;
	      maxsurf[icoord] -= Constants::twoPi;
	      midsurf[icoord] -= Constants::twoPi;
	    } else {
	      minsurf[icoord] += Constants::twoPi;
	      maxsurf[icoord] += Constants::twoPi;
	      midsurf[icoord] += Constants::twoPi;
	    }
	  }
	}
//
//  Update the global limits
//
	limits(0,icoord) = std::min(limits(0,icoord),minsurf[icoord]);
	limits(1,icoord) = std::max(limits(1,icoord),maxsurf[icoord]);
      }
      distlimits(0) = std::min(distlimits(0),mindist);
      distlimits(1) = std::max(distlimits(1),maxdist);
    }
  }
//
// put in a safety margin for the limits.
//
  distlimits(0) -= _perpdelta;
  distlimits(1) += _perpdelta;
  for(icoord=0;icoord<2;icoord++){
    limits(0,icoord) -= _surfdelta;
    limits(1,icoord) += _surfdelta;
    range(icoord) = _limits[1][icoord]-_limits[0][icoord];
    if(_wrapped[icoord]){
//
//  Force the wrapped coordinates to map evenly onto 2pi
//
      if(_range[icoord] >= Constants::twoPi*(_nbuckets[icoord]-1)/_nbuckets[icoord]){
//
//  If the limits are within 1 bucket of 2pi, stretch them to exactly 2pi
//
	limits(0,icoord) = 0.0;
	limits(1,icoord) = Constants::twoPi;
	range(icoord) = Constants::twoPi;
	nwrap(icoord) = _nbuckets[icoord];
      } else {
//
//  Otherwise, round up the range to give an integral number of buckets
//  from 0 to 2pi.
//
	double mfloat = Constants::twoPi*_nbuckets[icoord]/_range[icoord];
	nwrap(icoord) = int(mfloat);
	double newrange = Constants::twoPi*_nbuckets[icoord]/_nwrap[icoord];
	double epsi = (newrange - _range[icoord])/2.0;
	assert(epsi >= 0.0);
	range(icoord) = newrange;
	limits(0,icoord) -= epsi;
	limits(1,icoord) += epsi;
      }
    }
  }
}
//
//  Add a single element to the hash table
//
void
DetSurfaceSet::hashAddElem(DetElem* elem) const {
//
//  Get the physical outline
//
  std::vector<HepPoint> hpvec;
  elem->physicalOutline(hpvec);
  if(hpvec.size() == 0) return;
//
//  For now, just take the extrema of the outline in surface
//  coordinates.  An optimal algorithm would refine this
//  by testing each cell in the extrema to see if they actually
//  held part of the element.
//
  double mindist,maxdist;
  SurfacePoint minsurf,maxsurf,midsurf;
  elemSurfaceOutline(hpvec,minsurf,maxsurf,midsurf,mindist,maxdist);
//
//  Center the extrema in the limits for wrapped coordinates
//
  int icoord;
  for(icoord=0;icoord<2;icoord++)
    if(_wrapped[icoord])
      if(midsurf[icoord]<_limits[0][icoord]){
	minsurf[icoord] += Constants::twoPi;
	maxsurf[icoord] += Constants::twoPi;
	midsurf[icoord] += Constants::twoPi;
      } else if (midsurf[icoord]>_limits[1][icoord]){
	minsurf[icoord] -= Constants::twoPi;
	maxsurf[icoord] -= Constants::twoPi;
	midsurf[icoord] -= Constants::twoPi;
      }
//
//  Convert the extrema to hash table limits.
//
  int imin[2],imax[2];
  for(icoord=0;icoord<2;icoord++){
    imin[icoord] = hashIndex(minsurf[icoord],icoord);
    imax[icoord] = hashIndex(maxsurf[icoord],icoord);
  }
//
//  Loop over the cells, and insert
//
  for(int ix=imin[0];ix<=imax[0];ix++)
    for(int iy=imin[1];iy<=imax[1];iy++){
      int ihash = iHash(ix,iy);
      assert(ihash >= 0 && ihash < _nhash);
      _theHash[ihash].push_back(DetElemSurfRange(elem,mindist,maxdist));
    }
}
//
//  Make a list of elements from the cells covered by
//  a line segment defined by two surface points.
//
void
DetSurfaceSet::interElements(SurfacePoint* surfinter,
			     std::vector<const DetElem*>& intelements,
			     double* distrange) const{
//
//  Clear the list
//
  intelements.clear();
//
//  Convert the points to limits on the index ranges
//
  int ilim[2][2];
  int icoord;
  for(icoord=0;icoord<2;icoord++){
    if(_wrapped[icoord] && fabs(surfinter[0][icoord] -
			       surfinter[1][icoord])>Constants::pi){
//
//  Make sure wrapped points are on the same side of 2pi, and between
//  -2pi and 4pi (the implicit valid wrap region).
//
      if(surfinter[0][icoord]-surfinter[1][icoord]>Constants::pi){
	if(surfinter[1][icoord]<0)
	  surfinter[1][icoord] += Constants::twoPi;
	else
	  surfinter[0][icoord] -= Constants::twoPi;
      }
      else {
	if(surfinter[0][icoord] < 0)
	  surfinter[0][icoord] += Constants::twoPi;
	else
	  surfinter[1][icoord] -= Constants::twoPi;
      }
    }
    ilim[0][icoord] = hashIndex(std::min(surfinter[0][icoord],
                                         surfinter[1][icoord]),icoord);
    ilim[1][icoord] = hashIndex(std::max(surfinter[0][icoord],
                                         surfinter[1][icoord]),icoord);
  }
//
//  Now loop over the 2-D range, and build a list of elements
//
  for(int ix=ilim[0][0];ix<=ilim[1][0];ix++){
    for(int iy=ilim[0][1];iy<=ilim[1][1];iy++){
//
//  Convert this to a hash cell
//
      int ihash = iHash(ix,iy);
      if(ihash>=0){
//
//  Loop over the elements in this hash cell
//
	std::vector<DetElemSurfRange>& hashcell = _theHash[ihash];
	for ( std::vector<DetElemSurfRange>::const_iterator elemIter = hashcell.begin();
	      elemIter != hashcell.end();
	      ++elemIter ) {
	  const DetElemSurfRange& erange = *elemIter;
// test the depth; only insert elements which could be hit
	  if(distrange == 0 ||
	     !(distrange[1] < erange.minDist() ||
	       distrange[0] > erange.maxDist()) ){
//
//  Add the element to the list if it's not there already
//
            typedef std::vector<const DetElem*>::iterator iter_t;
            std::pair<iter_t,iter_t> p = 
                    std::equal_range(intelements.begin(),intelements.end(),erange.element());
	    if(p.first==p.second) {
	      intelements.insert(p.second,erange.element());
	    }
	  }
	}
      }
    }
  }
}

// integer hash, with final range check
int
DetSurfaceSet::iHash(int ixindex,int iyindex) const {
  int ix = hashBucket(ixindex,0);
  int iy = hashBucket(iyindex,1);
  if(ix >= 0 &&  iy >= 0 && ix < _nbuckets[0] && iy < _nbuckets[1])
    return iy + ix*_nbuckets[1];
  else
    return -1;
}
//
int
DetSurfaceSet::hashIndex(double val,int icoord) const{
  int ihash = int(floor(_nbuckets[icoord]*(val-_limits[0][icoord])/_range[icoord]));
  if(!_wrapped[icoord])
    ihash = std::min(std::max(ihash,0),_nbuckets[icoord]-1);
  return ihash;
}
//
//  Find the surface points corresponding to the intersection with the two surfaces
//  within the specified range (or the trajectory range, if the range is omitted).
//
bool
DetSurfaceSet::surfaceIntersections(const Trajectory* traj,
				    SurfacePoint* surfinter,
				    double* flightdist,
				    double* distrange,
				    double* range,
				    bool& between) const {
//
//  Test for the starting point to be between the cylinders.  This
//  also does a crude elimination of trajectorys completely
//  missing the surfaces.
//
  double startdist[2];
  double enddist[2];
  bool startbetween = betweenSurfaces(traj->position(range[0]),startdist);
  bool endbetween = betweenSurfaces(traj->position(range[1]),enddist);
  between = startbetween || endbetween;
// crude initial test; if the start isn't between and the distance to the
// closest surface is larger than the flight range, there's no hope of intersecting
// the surfaces.
  double pathlen = range[1] - range[0];
  if( (!startbetween) && (!endbetween) &&
      (std::min(fabs(startdist[0]),fabs(startdist[1])) > pathlen ||
       std::min(fabs(enddist[0]),fabs(enddist[1])) > pathlen) )
    return false;
//
//  Try forward intersecting the traj with the surfaces
//
  const Intersection inter[2] = { Intersection( *traj, *_surfaces[0] ),
				  Intersection( *traj, *_surfaces[1] ) };
  TrkErrCode iflag[2] = { inter[0].intersect(flightdist[0],trkOut,range),
			  inter[1].intersect(flightdist[1],trkOut,range) };
//
//  Now go through the cases: first, starting outside the surfaces
//
  bool returnval = false;
  int isurf=0, jsurf=0;
  if(!startbetween){
//  Both surfaces intersected: success
    if(iflag[0].success() && iflag[1].success())
      returnval = true;
    else if(iflag[0].success() || iflag[1].success()){
//  Only one intersected surface: this could be A) a looper or
//  B) a trajectory which ends between the surfaces.  Test these
      jsurf = (iflag[0].success() ? 1 : 0);
      isurf = (iflag[0].success() ? 0 : 1);
      if(endbetween){
//  Traj ends between the surfaces.  Set the 2nd point according to the
//  traj end.
	flightdist[jsurf] = range[1];
	returnval = true;
      } else {
//  Try re-intersecting the same surface a second time
	flightdist[jsurf] = flightdist[isurf]+_epsilon;
	iflag[jsurf] = inter[isurf].intersect(flightdist[jsurf],trkOut,range);
	if(iflag[jsurf].success())
	  returnval = true;
//  If re-intersecting didn't work, try a 'linear' intersection with the 2nd surface
	else {
	  double delta(0.0);
	  if(_surfaces[jsurf]->
	     distTo(traj->position(flightdist[isurf]),
		    traj->direction(flightdist[isurf]),
		    delta,
		    DetSurface::closest) != DetSurface::nointersect){
	    flightdist[jsurf] = flightdist[isurf] + delta;
	    returnval = true;
// error; do the best we can with the point we have
	  } else {
	    ErrMsg(routine) << "DetSurfaceSet: inconsistent surface intersections" << endmsg;
	    flightdist[jsurf] = range[jsurf];
	    returnval = true;
	  }
	}
      }
    }
  } else {
//  Starting between the surfaces: this establishes one point.  Check the
//  intersections for the remaining point.
    if(iflag[0].success() && iflag[1].success()){
//  Both surfaces intersected: take the second point as the one with the
//  shortest flight distance
      jsurf = (flightdist[0]>flightdist[1] ? 0 : 1);
      flightdist[jsurf] = range[0];
      returnval = true;
    } else if(iflag[0].success() || iflag[1].success()){
//  Only one intersected surface: this defines the other point
      jsurf = (iflag[0].success() ? 1 : 0);
      flightdist[jsurf] = range[0];
      returnval = true;
//  No intersection: maybe the end point is between the surfaces too.  If so,
//  it defines the second point.
    } else {
      if(!endbetween)
	ErrMsg(routine) << "DetSurfaceSet: inconsistent surface intersections" << endmsg;
      flightdist[0] = range[0];
      flightdist[1] = range[1];
      returnval = true;
    }
  }
// surface points use the reference surface; this also computes the perp distance
  Hep3Vector norm;
  for(isurf=0;isurf<2;isurf++)
    distrange[isurf] = -project(traj->position(flightdist[isurf]),norm,surfinter[isurf]);
// check for 'reversed' trajectories
  if(distrange[0]>distrange[1]){
    const double lowdist = distrange[1];
    distrange[1] = distrange[0];
    distrange[0] = lowdist;
  }
  return returnval;
}
//
//  Test for a point to be between the surfaces
//
bool
DetSurfaceSet::betweenSurfaces(const HepPoint& point,double* dist) const {
//
//  Find the normals from the confining surfaces to this point.
//  Note the FUNKY SIGN CONVENTION relative to segminmax
//
  Hep3Vector normvec;
//   for(int isurf=0;isurf<2;isurf++)
//     dist[isurf] = -_surfaces[isurf]->normTo(point,normvec);
  dist[0] = -_surfaces[0]->normTo(point,normvec);
  dist[1] = -_surfaces[1]->normTo(point,normvec);
  return dist[0]*dist[1] < 0.0;
}
//
//  printout and diagnostics
//
void
DetSurfaceSet::printAll(ostream& os) const {
  if(isReady()) {
    os << "Detector Surface Set " << _dsname << " " << _dsnum
       << " is bounded by the following surfaces: " << endl;
    os << "Inner surface = ";
    _surfaces[0]->printAll(os);
    os << "Outer surface = ";
    _surfaces[1]->printAll(os);
    os << "Hash table has " << _nbuckets[0] << " buckets in view 0 and "
       << _nbuckets[1] << " buckets in view 1" << endl;
    os << "Surface coordinate view 0 limits are from "
       << _limits[0][0] << " to " << _limits[1][0];
    if(_wrapped[0])
      os << " and is wrapped " << endl;
    else
      os << " and is NOT wrapped " << endl;
    os << "Surface coordinate view 1 limits are from "
       << _limits[0][1] << " to " << _limits[1][1];
    if(_wrapped[1])
      os << " and is wrapped " << endl;
    else
      os << " and is NOT wrapped " << endl;  
// Build up a list of all the elements
    DetElemList allElemList;
    listAllElements(allElemList);
    os << "There are " << allElemList.size() 
       << " elements controled by this set, assigned to hash bins as follows "
       << endl;
// go through the list and print the hash table entries for each element
    DetElemList::const_iterator iter = allElemList.begin();
    DetElem* elem = 0;
    while ( iter != allElemList.end() ) {
      elem = *iter++;
      printElement(elem,os);
    }
  } else
    os << "Detector Surface Set " << _dsname << " " << _dsnum
       << " is NOT READY for intersections" << endl;
}

// print the hash cells which contain a particular element
void
DetSurfaceSet::printElement(const DetElem* elem,ostream& os) const {
  elem->print(os);
  os << "Contained in the following cells: ";
//  os << setw(3) << setprecision(3);
  for(int ihash=0;ihash<_nhash;ihash++){
    unsigned nelem = _theHash[ihash].size();
    for(unsigned ielem=0;ielem<nelem;ielem++){
      std::vector<DetElemSurfRange>& aHashBucket= _theHash[ihash];
      if(aHashBucket[ielem].element() == elem) {
	int ix = ihash/_nbuckets[1];
	int iy = ihash - ix*_nbuckets[1];
	os << ":" << ix << "," << iy;
	break;
      }
    }
  }
  os << endl;
}

double
DetSurfaceSet::project(const HepPoint& point,Hep3Vector& norm,
		       SurfacePoint& surfpoint) const {
  return _refsurface->normalTo(point,norm,surfpoint);
}
