//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetVolumeElem.cc,v 1.27 2004/08/06 05:58:33 bartoldu Exp $
//
// Description:
//	Class DetVolumeElem. See header for information.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Stephen J. Gowdy	        Originator
//	Gautier Hamel de Monchenault	Originator
//
// Copyright Information:
//	Copyright (C) 1997	University of Edinburgh
//	Copyright (C) 1997	CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "BaBar/Constants.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DetectorModel/DetVolumeElem.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------
#include <assert.h>
#include <vector>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Geometry/Transformation.h"
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetSurface.hh"
#include "DetectorModel/DetVolSideIntersection.hh"
#include "DetectorModel/DetVolumeType.hh"
#include "DetectorModel/Intersection.hh"
using std::endl;
using std::ostream;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DetVolumeElem::DetVolumeElem( DetVolumeType* itsType, const char* name, int id,
			      const HepTransformation& theAlignment )
  : DetElem( itsType, name, id ),
    _sides( new std::vector< DetSurface* > )
{
  myTransf() = new HepTransformation( theAlignment );
  createCache();
}

//--------------
// Destructor --
//--------------
DetVolumeElem::~DetVolumeElem()
{
  delete myTransf();

  std::vector< DetSurface* >::iterator siter = _sides->begin();
  while( siter != _sides->end() ) {
    DetSurface* aSurf = *siter++;
    if( 0 != aSurf ) {
      delete aSurf;
      aSurf = 0;
    }
  }
  _sides->clear();
  delete _sides;
  _sides = 0;
}

//-------------
// Methods   --
//-------------
    
//-------------
// Operators --
//-------------
    
//-------------
// Selectors --
//-------------
int 
DetVolumeElem::intersect( const Trajectory* traj,
			  DetIntersection& dinter ) const
{
  const double minCurvature=0.0001*0.25*Constants::pi;

  // converted from RWTValSortedVector
  std::vector< DetVolSideIntersection > ilist;

  //search for sides up to 1/2 turn away
  double range[2];
  range[0]=dinter.pathrange[0];
  range[1]=dinter.pathrange[1];
  double k=fabs(traj->curvature(dinter.pathrange[0]));
  double vrho = (traj->direction(dinter.pathrange[0])).perp();
  if(k<minCurvature*vrho) {
    range[0]-=0.25*Constants::pi/minCurvature;
    range[1]+=0.25*Constants::pi/minCurvature;
  }else {
    range[0]-=0.25*Constants::pi*vrho/k;
    range[1]+=0.25*Constants::pi*vrho/k;
  }

  //find intersections with sides
  sideIntersect( traj, ilist, dinter.pathlen, range);

  //classify intersections 
  if(ilist.size()==2) {

    double entrance=ilist[0].pathLength();
    double exit=ilist[1].pathLength();
    
    dinter.flag[0]=-1301;	//doesn't get to element 
    if(entrance>dinter.pathrange[1]) return false;
    
    dinter.flag[0]=-1303;	// completely past element
    if(exit<dinter.pathrange[0]) return false;
    
   if(entrance<dinter.pathrange[0]
       &&
       exit>dinter.pathrange[1]) {
      entrance=dinter.pathrange[0];
      exit=dinter.pathrange[1];
      double path=0.5*(entrance+exit);
      dinter=DetIntersection
      	((DetElem*)this,traj,path,entrance,exit);
     
       dinter.flag[0]=-1;	// completely in element
       dinter.flag[1]=-1;
       return true;
    }
    
    
    if(entrance>=dinter.pathrange[0]
       && 
       exit<=dinter.pathrange[1]) {
      
      
      //thru going trajectory
      
      double path     = 0.5*(entrance+exit);
      dinter = DetIntersection( (DetElem*)this, 
				traj, path, entrance, exit);
      dinter.flag[0]  = ilist[0].side();
      dinter.flag[1]  = ilist[1].side();
      
    }else {
      
      if(exit>dinter.pathrange[1]) {
	
	//doesn't get all the way through element
	exit=dinter.pathrange[1];
	double path=0.5*(entrance+exit);
	dinter=DetIntersection
	  ((DetElem*)this,traj,path,entrance,exit);
	dinter.flag[0]=ilist[0].side();
	dinter.flag[1]=-1;
	
      }
      else{
	
	//trajectory starts in element
	
	entrance=dinter.pathrange[0];
	double path=0.5*(entrance+exit);
	dinter=DetIntersection
	  ((DetElem*)this,traj,path,entrance,exit);
	dinter.flag[0]=-1;
	dinter.flag[1]=ilist[1].side();
	
      }
	
    }
    return true;
  }

  //failed to find both sides or found too many sides
  dinter.flag[0]=-1304;	
  dinter.flag[1]=ilist.size();
  return false;
}
    
void
DetVolumeElem::sideIntersect
  ( const Trajectory* traj,
    std::vector< DetVolSideIntersection >& ilist,
    double distance,
    double* range ) const
{
  int jflag = 0;
  const double epsilon = 1.0e-05;
  DetVolumeType* volType = (DetVolumeType*) _dtype;
  SurfacePoint surfcoord;
  size_t side=0, nSides = _sides->size();
  for( side=0; side<nSides; side++)
    {
      //
      //  Use the initial values to define the starting point and the 
      //  search range
      //
      double flightdist = distance;
      double flightrange[2];
      flightrange[0] = range[0];
      flightrange[1] = range[1];
      Intersection intersection( *traj, *(*_sides)[side] );
      TrkErrCode iflag = intersection.intersect( flightdist, 
						 surfcoord,
						 trkOut,
						 flightrange );
      while(iflag.success())
	{
	  TwoDCoord detcoord(surfcoord.array());
	  if( !(volType->insideLimitsOf( side, surfcoord )) )
	    {
	      // not inside the physical boundaries for this plane,
	      // go to the next intersection
	      flightrange[0] = flightdist + epsilon;
	      iflag =  intersection.intersect( flightdist, 
					       surfcoord,
					       trkOut,
					       flightrange );
	    } 
	  else
	    {
	      ilist.push_back(DetVolSideIntersection(flightdist,side,detcoord));
	      break;
	    }
	}
    }
  std::sort(ilist.begin(), ilist.end() );
}
    
void
DetVolumeElem::physicalOutline(std::vector<HepPoint>& pvec) const {
  pvec.clear();
 
  const std::vector< TypeCoord* >* theOutline
    = detectorType()->pointStore();
  //
  //  Loop over the corners
  //
  int npoints = theOutline->size();
  for( int icorn=0; icorn<npoints; icorn++ )
    {
      //
      //  Turn them into space points
      //
      ThreeDCoord* typeCoord = (ThreeDCoord*)(*theOutline)[icorn];
      HepPoint typePoint( (*typeCoord)[0],
			  (*typeCoord)[1],
			  (*typeCoord)[2] );
      pvec.push_back( _etrans->transFrom( typePoint ) );
    }
}


void 
DetVolumeElem::print(ostream& o) const 
{
  DetElem::print( o );
}

void 
DetVolumeElem::printAll(ostream& o) const 
{
  DetElem::print( o );
  o << "Reference Surface = DetPlane" << endl;
  o << " Center point " <<_etrans->origin() <<  endl;
  o << " x axis       " <<_etrans->unit(I_x) <<  endl;
  o << " y axis       " <<_etrans->unit(I_y) <<  endl;
  o << " z axis       " <<_etrans->unit(I_z) <<  endl;
  size_t side=0, nSides = _sides->size();
  for ( side=0; side<nSides; side++ )
    {
      o << "Side " << side << endl;
      o << (*_sides)[side]->centerPoint() << endl;
      o << (*_sides)[side]->basis(0) << endl;
      o << (*_sides)[side]->basis(1) << endl;
      o << (*_sides)[side]->basis(2) << endl;
    }
}

//-------------
// Modifiers --
//-------------
void
DetVolumeElem::createCache()
{
  DetVolumeType* volType = (DetVolumeType*) detectorType();

  size_t side=0, nSides = volType->_sides->size();
  for( side = 0; side<nSides; side++ )
    {
      HepTransformation tran( transf() );
      DetSurface* thisSurface = (*volType->_sides)[side];
      _sides->push_back( thisSurface->copyOf() );
      assert( _sides->size() == side+1 );
      tran *= *( thisSurface->transform() );
      *( (*_sides)[side]->transform() ) = tran;
    }
}

void
DetVolumeElem::updateCache()
{
  DetVolumeType* volType = (DetVolumeType*) detectorType();

  size_t side=0, nSides = _sides->size();
  for( side = 0; side<nSides; side++ )
    {
      HepTransformation tran( transf() );
      DetSurface* thisSurface = (*volType->_sides)[side];
      tran *= *( thisSurface->transform() );
      *( (*_sides)[side]->transform() ) = tran;
    }
}

//		-----------------------------------------------
// 		-- Static Data & Function Member Definitions --
//		-----------------------------------------------
ostream&
operator<<( ostream& o, const DetVolumeElem& a )
{
  a.print(o); 
  return o;
}

//		-------------------------------------------
// 		-- Protected Function Member Definitions --
//		-------------------------------------------
HepPoint
DetVolumeElem::coordToPoint( const TypeCoord* aCoord ) const
{
  HepPoint aPoint( (*aCoord)[0], (*aCoord)[1], (*aCoord)[2] );
  return transform().transFrom( aPoint );
}

//		-----------------------------------------
// 		-- Private Function Member Definitions --
//		-----------------------------------------

//		-----------------------------------
// 		-- Internal Function Definitions --
//		-----------------------------------
