//
//  Code for the DetSurfaceType class
//
#include "BaBar/BaBar.hh"
#include "ErrLogger/ErrLog.hh"
#include <math.h>
#include "CLHEP/Utilities/CLHEP.h"
#include "DetectorModel/DetSurfaceType.hh"
#include "DetectorModel/DetSurface.hh"
#include <vector>
using std::vector;
//
DetSurfaceType::DetSurfaceType() :
  DetType("Unknown",-1)
{;}
//
DetSurfaceType::DetSurfaceType(const char* name,int ityp) :
  DetType(name,ityp)
{;}

DetSurfaceType::DetSurfaceType(const char* name,int ityp,
			       std::vector< TwoDCoord >& theCorners ) :
  DetType(name,ityp)
{
  // Add the points to the pointStore list from the
  // DetType class
  std::vector< TypeCoord* >* thePointStore = myPointStore();

  size_t i=0, nCorners = theCorners.size();
  for ( i=0; i<nCorners; i++ )
    {
      TypeCoord* thePoint = &theCorners[i];
      thePointStore->push_back( thePoint->copyOf() );
      // It also makes up the outline list used to draw it.
      myOutline()->push_back( (*thePointStore)[ i ] );
    }
  myOutline()->push_back( (*myPointStore())[ 0 ] );
}

DetSurfaceType::~DetSurfaceType()
{
}

//
//  Symmetric implementation: classes which care about
//  the coordinate direction cosines separately should override
//  this function.
//
double
DetSurfaceType::effectiveThickness(const TwoDCoord* point,
				double dircos1,
				double dircos2) const {
  return thickness(point)/sqrt(1.0 - (sqr(dircos1)+sqr(dircos2)));
}

bool
DetSurfaceType::insideLimitsOf( const TwoDCoord& thisPoint ) const
{

  const std::vector< TypeCoord* >* thePointStore = pointStore();
  unsigned numberOfCorners = thePointStore->size();
  if(numberOfCorners == 0) return false;
  unsigned i=0;

  // Loop over all the edges of this side and see if the point is inside
  //   (inside means the same side as the origin)
  const TwoDCoord* point1;
  const TwoDCoord* point2;
  for ( i=0; i < ( numberOfCorners - 1 ); i++ )
    {
      point1 = (const TwoDCoord*) (*thePointStore)[i];
      point2 = (const TwoDCoord*) (*thePointStore)[i+1];
      if ( ! insideLine( thisPoint, *point1,
			 *point2 ) ) return false;
    }
  // To check the last edge has a special case
  point1 = (const TwoDCoord*) (*thePointStore)[numberOfCorners-1];
  point2 = (const TwoDCoord*) (*thePointStore)[0];
  if ( ! insideLine( thisPoint, *point1,
		     *point2 ) ) return false;

  return true;
}

bool
DetSurfaceType::insideLine( const TwoDCoord& toTest,
			    const TwoDCoord& p1,
			    const TwoDCoord& p2 ) const
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

  return answer;
}
