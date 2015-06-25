//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetSimpleVolumeType.cc,v 1.20 2003/01/27 22:31:36 ryd Exp $
//
// Description:
//	Class DetSimpleVolumeType
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
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/ErrLog.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "BTrk/DetectorModel/DetSimpleVolumeType.hh"

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
using std::vector;

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BTrk/BbrGeom/Transformation.h"
#include "BTrk/DetectorModel/DetPlane.hh"
#include "BTrk/DetectorModel/DetSurface.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructor  --
//----------------
DetSimpleVolumeType::DetSimpleVolumeType( const char* typeName, int typeId,
					  const DetMaterial* theMaterial,
					  std::vector< HepPoint >& theCorners )
 : DetVolumeType( typeName, typeId ), _material( theMaterial )
{
  // Add the points to the pointStore list from the
  // DetType class
  std::vector< TypeCoord* >* thePointStore = myPointStore();

  size_t i=0, nCorners = theCorners.size();
  for ( i=0; i<nCorners; i++ )
    {
      HepPoint* thePoint = &theCorners[i];
      thePointStore->push_back( new ThreeDCoord( thePoint->x(),
					      thePoint->y(),
					      thePoint->z() ) );
    }

  size_t side=0, corner=0, nEdges=0, nSides=0, nEdgesThisSide=0;

  // nEdges is found by dividing the number of corners by two
  nEdges = nCorners/2;
  nSides = nEdges + 2;
  assert( nCorners == nEdges*2 );
  std::vector< HepPoint* > sideCorners( nEdges );

  // Reshape the sideCorners defined in parent class to have correct
  // number of lists
  mySideCorners()->resize( nEdges + 2 );

  // Make sides from corners..
  // It also makes up the outline list used to draw it.
  // It will draw the two main sides and then lines between
  // each corner of these sides.
  for ( side=0; side<nSides; side++ )
    {
      switch(side)
	{
	case near:
	  for ( corner=0; corner<nEdges; corner++ )
	    {
	      // Order corners backwards to get normal pointing inwards
	      sideCorners[nEdges-corner-1] = &theCorners[corner];
	      myOutline()->push_back( (*myPointStore())[corner] );
	    }
	  myOutline()->push_back( (*myPointStore())[ 0 ] );
	  myOutline()->push_back( (ThreeDCoord*)0 );
	  nEdgesThisSide = nEdges;
	  break;

	case far:
	  for ( corner=nEdges; corner<(nEdges*2); corner++ )
	    {
	      sideCorners[corner-nEdges] = &theCorners[corner];
	      myOutline()->push_back( (*myPointStore())[ corner ] );
	    }
	  myOutline()->push_back( (*myPointStore())[ nEdges ] );
	  myOutline()->push_back( (ThreeDCoord*)0 );
	  nEdgesThisSide = nEdges;
	  break;

	default:
	  if ( side == nEdges+1 )
	    {
	      // Special case of last side
	      sideCorners[0] = &theCorners[nEdges-1];
	      sideCorners[1] = &theCorners[0];
	      sideCorners[2] = &theCorners[nEdges];
	      sideCorners[3] = &theCorners[nEdges*2-1];
	    }
	  else
	    {
	      sideCorners[0] = &theCorners[side-2];
	      sideCorners[1] = &theCorners[side-1];
	      sideCorners[2] = &theCorners[side+nEdges-1];
	      sideCorners[3] = &theCorners[side+nEdges-2];
	    }
	  myOutline()->push_back( (*myPointStore())[ side-2 ] );
	  myOutline()->push_back( (*myPointStore())[ side+nEdges-2 ] );
	  myOutline()->push_back( (ThreeDCoord*)0 );
	  nEdgesThisSide = 4;
	  break;
	}
      // Resize the container of side corners for this side
      (*mySideCorners())[side].resize( nEdgesThisSide );

      Hep3Vector sideOne = *sideCorners[0] - *sideCorners[1],
	sideTwo = *sideCorners[2] - *sideCorners[1],
	normal = sideOne.cross(sideTwo);

      Hep3Vector centre;

      for ( i=0; i<nEdgesThisSide; i++ )
	centre += Hep3Vector( sideCorners[i]->x(),
			      sideCorners[i]->y(),
			      sideCorners[i]->z() );

      centre *= 1./nEdgesThisSide;

      std::vector< DetSurface* >* theSides = mySides();
      theSides->push_back( new DetPlane( HepTransformation( centre, normal ) ) );

      for ( i=0; i<nEdgesThisSide; i++ )
	{
	  // Fill up container will all the local side corners
	  SurfacePoint* aPoint = new SurfacePoint;
	  Hep3Vector dummy;
	  (*mySides())[side]->normalTo( *sideCorners[i],
					dummy,
					*aPoint );
	  (*mySideCorners())[side][i] = aPoint;
	}
    }
}



//--------------
// Destructor --
//--------------
DetSimpleVolumeType::~DetSimpleVolumeType()
{
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
    
//-------------
// Modifiers --
//-------------

//		-----------------------------------------------
// 		-- Static Data & Function Member Definitions --
//		-----------------------------------------------

//		-------------------------------------------
// 		-- Protected Function Member Definitions --
//		-------------------------------------------
bool
DetSimpleVolumeType::insideLine( const SurfacePoint& toTest,
				 const SurfacePoint& p1,
				 const SurfacePoint& p2 ) const
{
  bool answer=false;
  double m=0, c=0;
  double mt=0;
  double xpt=0, ypt=0;

  // See if p1 and p2 is a vertical line
  if ( p2[0]==p1[0] )
    {
      // Test to see if toTest is nearer origin (0,0) than line
      if ( ( p1[0] - toTest[0] ) * p1[0] >= 0 ) answer = true;
    }
  else
    {
      // Find line between the two reference points
      m = ( p2[1]-p1[1] ) / (p2[0]-p1[0]);
      c = p1[1] - m*p1[0];

      // Find the line between the test point and the origin
      mt = toTest[1]/toTest[0];

      if( fabs(mt-m) == 0. )
	{
	  answer = true;
	}
      else
	{

	  // Find where lines cross
	  xpt = c/(mt-m);

	  // If this point is further from the origin then
	  // the point is inside the line (unless the point is
	  // on the other side of the origin, then it is defined
	  // as inside the line.)
	  if ( ( xpt >= 0 && toTest[0] >= 0 && xpt >= toTest[0] ) ||
	       ( xpt < 0 && toTest[0] < 0 && xpt < toTest[0] ) ||
	       ( xpt < 0 && toTest[0] >= 0 ) ||
	       ( xpt >= 0 && toTest[0] < 0 ) ) answer=true;
	}
    }

  return answer;
}

//		-----------------------------------------
// 		-- Private Function Member Definitions --
//		-----------------------------------------

//		-----------------------------------
// 		-- Internal Function Definitions --
//		-----------------------------------
