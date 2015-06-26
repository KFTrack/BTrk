// ---------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetSurfaceSet.hh,v 1.24 2004/12/14 07:10:18 bartoldu Exp $
//
//  Description:
//  Special form of a DetSet where the elements are arranged approximately
//  along a surface.  Any type of surface (DetSurface subclass) can be used.
//  The surface arrangement of the elements is exploited in building a hash
//  table of their position on the surface, allowing faster navigation.  The
//  user must supply a seed surface and functions for creating a family of
//  similair surfaces from the seed as part of the construtor.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/1/97
//----------------------------------------------------------------------------
#ifndef DETSURFACESET_HH
#define DETSURFACESET_HH
//
//  Includes
//
#include "BTrk/DetectorModel/DetSet.hh"
#include "BTrk/DetectorModel/DetSurface.hh"
#include "BTrk/DetectorModel/DetElemSurfRange.hh"
#include <vector>

class HepPoint;
//
//  Function for list searches
//
bool findElemByPtr(DetElem*,void*);
//
//  Define the class
//
class DetSurfaceSet : public DetSet {
public:
  enum bounds { inner=0,outer=1};
//  Unique constructor; this needs pointers to 2 functions and a reference surface
  DetSurfaceSet(const char* name,int IDNum, // usual DetSet stuff
		const DetSurface& seed, // reference surface
		int nxash=20,int nyhash=20); // by default, 20X20 hash table
//  Destructor
  virtual ~DetSurfaceSet();
//
//  Overwrite default functions with the improved functions.  See the base class
//  DetSet for a description of these
//
  virtual bool firstIntersection(const Trajectory* traj,DetIntersection& next,
				 double* myrange=0) const;
//
  virtual void intersection(std::vector<DetIntersection>&,
			    const Trajectory*,double* myrange=0,
			    bool clear=true) const;

//
//  Overwrite the alignment functions, to keep the DetSurfaces OK.
//
  void applyGlobal(const DetAlignElem&);// apply global alignment
  void removeGlobal(const DetAlignElem&);// unapply global alignment
// overwrite printAll, to dump the hash table
  void printAll(std::ostream& os ) const;
// special diagnostic printout
  void printElement(const DetElem* elem,std::ostream& os) const;
//
  bool isReady() const;
  void makeReady() const;
// access to the surfaces
  const DetSurface* referenceSurface() const { return _refsurface; }
  const DetSurface* boundingSurface(bounds isurf) const { return _surfaces[isurf]; }
protected:
// allow overriding the surface coordinates used for the hashing.  This
// function should transform the point in space to a surface point and
// the perpendicular distance of that point from the REFERENCE SURFACE.
// Default implementation is to use the normalTo function of the reference
// surface.
  virtual double project(const HepPoint&,CLHEP::Hep3Vector&,SurfacePoint&) const;
private:
//
//  Data members
//
  DetSurface* _refsurface; //copy of the reference surface
  DetSurface* _surfaces[2]; // Surfaces circumscribing all elements
  int _nbuckets[2]; // number of buckets in the hash table in each view
  int _nhash; // total number of hash cells
  double _limits[2][2]; // define the surface coordinate range of the set
  double _range[2]; // range implied by the above limits
  double _distlimits[2]; // define the out-of-surface limits of the set
  bool _wrapped[2]; // store whether coordinates are wrapped
  int _nwrap[2];  // for wrapped coordinates, the number of buckets from 0->2pi

  mutable std::vector< std::vector<DetElemSurfRange> > _theHash; 
                                           // hash table of element lists,
                                           // PtrVector of owned pointers

//
//  Private functions
//
  bool withinRange(const SurfacePoint&) const;//range test function
//  find elements in cells between 2 surface points
  void interElements(SurfacePoint*,
		     std::vector<const DetElem*>&,
		     double* distrange=0) const;
//  point is between the surfaces
  bool betweenSurfaces(const HepPoint&,double*) const;
//
//  Functions to cast-off const for cache items
//
  DetSurface*& surfaces(int isurf) const {
    return (DetSurface*&)_surfaces[isurf]; }
  DetSurface*& refSurface() const {
    return (DetSurface*&)_refsurface; }
  double& limits(int iminmax,int icoord) const {
    return (double&)_limits[iminmax][icoord]; }
  double& distlimits(int iminmax) const {
    return (double&)_distlimits[iminmax]; }
  double& range(int icoord) const {
    return (double&)_range[icoord];}
  int& nwrap(int icoord) const {
    return (int&)_nwrap[icoord];}

//
//  Functions for building the hash table
//
  void buildHash() const; // function to build hash table
  void globalLimits(DetElemList&) const; // global limits
  void hashAddElem(DetElem*) const; // add a single element
//
//  Function to translate an elements perimeter (physicalOutline) to 
//  surface coordinates limits.  This also computes the perp. distance
//  extrema.
//
  void elemSurfaceOutline(std::vector<HepPoint>&,SurfacePoint&,
			  SurfacePoint&,SurfacePoint&,
			  double&,double&) const;
//
//  Function to intersect a trajectory with the surfaces
//
  bool surfaceIntersections(const Trajectory*,SurfacePoint*,
			    double*,double*,double* frange,bool&) const;
//
//  Hashing helper functions, inlined for speed
//
  int hashIndex(double val,int icoord) const;
  int iHash(int,int) const;
  int hashBucket(int index,int icoord) const { // assure the index limits
    if(_wrapped[icoord]){ // try to push wrapped coordinates out-of-range back in
      if(index < 0)
	index += _nwrap[icoord];
      else if(index >= _nbuckets[icoord])
	index -= _nwrap[icoord];
    }
    if(index >= 0 && index < _nbuckets[icoord])
      return index;
    else
      return -1;
  }
// prohibit
  DetSurfaceSet& operator = (const DetSurfaceSet&);
  DetSurfaceSet(const DetSurfaceSet&);
};
#endif
