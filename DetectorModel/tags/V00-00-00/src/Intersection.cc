// $Id: Intersection.cc,v 1.34 2004/10/11 19:08:05 brownd Exp $
//
// Coded stolen from gismo 2 
//    and modified for BaBar by Gautier Hamel de Monchenault
#include "BaBar/BaBar.hh"
#include "DetectorModel/Intersection.hh"
#include <math.h>
#include "BaBar/Constants.hh"
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
using std::cerr;
using std::cout;
using std::endl;

static const double minstep = 1e-5;       // step size to exit iterations
static const double maxrange = 5000.; // maximum range for track intersections
static const unsigned maxit = 100; // maximum number of steps over the range
static const double intolerance = 0.01;
static const double invsqrt2 = 0.7071067812;

TrkErrCode 
Intersection::intersect( double& fdist,
			 trkDirection tdir,double* frange) const {
//
//  First, establish the valid limits for the search; these are the combination of
//  the explicit input range (if provided) and the trajectory limits.
//
  double f1 = tdir == trkOut ? std::max(fdist,-maxrange) : std::min(fdist,maxrange);
  double fmax(maxrange);
  double fmin(-maxrange);
  if(frange){
    f1 = tdir == trkOut ? std::max(f1,frange[0]) : std::min(f1,frange[1]);
    fmax = std::min(fmax,frange[1]);
    fmin = std::max(fmin,frange[0]);
  } else {
    f1 = tdir == trkOut ? std::max(f1,traj.lowRange()) : std::min(f1,traj.hiRange());
    fmax = std::min(fmax,traj.hiRange());
    fmin = std::max(fmin,traj.lowRange());
  }
// get an approximate intersection position that's forward of the
// starting point
  DetSurface::intertype inter = findProximity(f1,fmin,fmax,tdir);
// iterate to convergence
  unsigned iter(0);
  double delta = f1-fmin;
  TrkErrCode retval(TrkErrCode::succeed);
  HepPoint   p;
  Hep3Vector v;
  while(inter != DetSurface::nointersect) {
// update the position and limits
    traj.getInfo( f1, p, v );
// here we're looking for the closest solution
    inter = surf.distTo(p,v,delta,DetSurface::closest);
// update
    f1 += delta;
// convergence test
    if(fabs(delta) <= minstep )break;
// count the number of iterations
    iter++;
    if(iter>maxit){
      retval = TrkErrCode(TrkErrCode::fail,59,"Step number limit reached");
      break;
    }
  }
// convergence: test to see what happened
  if(inter == DetSurface::intersect){
// last check on the range
    if(f1>=fmin && f1 <= fmax)
      fdist = f1; // done
    else
      retval = TrkErrCode(TrkErrCode::fail,55,"No intersection within range");
  } else if(inter == DetSurface::localmin)
    retval = TrkErrCode(TrkErrCode::fail,58,"Local Minimum");
  else
    retval = TrkErrCode(TrkErrCode::fail,59,"No Intersection");
#ifdef VERIFYINTERSECTION
// call the old code to test the new
  double oldf(f0);
  SurfacePoint olduv;
  TrkErrCode oldinter = oldintersect(oldf,olduv,tdir,frange);
  if(retval.failure() && oldinter.success()){
// verify it
    p = traj.position(oldf);
    if(surf.surfacePoint(p,olduv,minstep)==0 &&
       oldf >= fmin && oldf <= fmax){
      cerr << "Error: valid old intersection found without new" << endl;
    } else 
       cerr << "Error: invalid old intersection found" << endl;
  } else if(retval.success() && oldinter.failure()){
// test it again
    p = traj.position(fdist);
    Hep3Vector norm;
    double sdist = surf.normTo(p,norm);
    if(fdist>=fmin && fdist <= fmax && fabs(sdist)<minstep){
      cout << "Found new intersection missed by old code" << endl;
    } else
      cerr << "Error: found bad new intersection, distance = "
	   << sdist << endl;
  } else if(retval.success() && oldinter.success() &&
	    fabs(fdist-oldf)>minstep){
// test them both
    p = traj.position(fdist);
    Hep3Vector norm;
    double newdist = surf.normTo(p,norm);
    p = traj.position(oldf);
    double olddist = surf.normTo(p,norm);
    if(fdist>=fmin && fdist <= fmax && fabs(newdist)<minstep &&
       oldf>=fmin && oldf <= fmax && fabs(olddist)<minstep){
      if(fdist<oldf)
	cout << "Found new intersection before old" << endl;
      else {
	 cerr << "Error: Found old intersection "
	      << oldf << " before new " << fdist
	      << " Surface,Trajectory = " << endl;
	 surf.printAll();
	 traj.printAll();
      }
    } else if (fdist>=fmin && fdist <= fmax && fabs(newdist)<minstep)
      cerr << "Found bad old intersection" << endl;
    else
      cerr << "Error: Found bad new intersection" << endl;
  }
#endif
  return retval;
}

TrkErrCode 
Intersection::intersect( double& fdist, SurfacePoint& uv,
			 trkDirection tdir,double* frange) const {
  TrkErrCode retval = intersect(fdist,tdir,frange);
// compute the surface point
  if(retval.success()){
    HepPoint p = traj.position(fdist);
    if(surf.surfacePoint(p,uv,minstep)!=0)
      retval = TrkErrCode(TrkErrCode::fail,57,"Point not on surface");
  }
  return retval;
}

TrkErrCode 
Intersection::oldintersect( double& f, SurfacePoint& uv,
			 trkDirection tdir,double* frange) const
{
  Hep3Vector n;
  HepPoint   p;
  Hep3Vector v;
//
//  First, establish the valid limits for the search; these are the combination of
//  the explicit input range (if provided) and the trajectory limits.
//
  double f1 = tdir == trkOut ? std::max(f,-maxrange) : std::min(f,maxrange);
  double fmax(maxrange);
  double fmin(-maxrange);
  if(frange){
    f1 = tdir == trkOut ? std::max(f1,frange[0]) : std::min(f1,frange[1]);
    fmax = std::min(fmax,frange[1]);
    fmin = std::max(fmin,frange[0]);
  } else {
    f1 = tdir == trkOut ? std::max(f1,traj.lowRange()) : std::min(f1,traj.hiRange());
    fmax = std::min(fmax,traj.hiRange());
    fmin = std::max(fmin,traj.lowRange());
  }
  if( (tdir == trkOut ? f1 >= fmax : f1 <= fmin) ) return 
					    TrkErrCode(TrkErrCode::fail,2);
  f = f1;
  traj.getInfo( f, p, v );
//
//  Perform a coarse step search, to find points which stradle the surface.
//  First determine the step size.  This combines the information of
//  how far we are from the surface with how fast the trajectory and
//  the surface are changing.
//
  double tdist = tdir == trkOut ? fabs(fmax-f1) : fabs(f1-fmin);
  double d1(surf.normalTo( p, n, uv ));
  double scurve = surf.curvature(uv);
  double sstep = (scurve > 0) ? 0.1*Constants::pi/scurve : tdist;
  double tstep = traj.distTo1stError(f,intolerance,tdir);
  double delta = std::max(std::min(fabs(sstep),fabs(tstep)),tdist/maxit);
//
//  Arrange the step so that the full range is an integral number of steps
//
  int niter = std::max(1,(int)rint(tdist/delta));
  delta = tdir == trkOut ? (fmax-f1)/niter : (fmin-f1)/niter ;
//
//  Loop till we converge or exceed the valid range or run out of patience
//
  double d2(d1);
  double f2(f1);
  int iter=0;

  while( iter < niter ) {
    f2 = f1 + delta;
    p = traj.position( f2 );  
    d2 = surf.normalTo( p, n, uv );
    if(d1*d2 < 0 ) {
//
//  We now have points on either side of the surface;
//  Find a precise solution through interpolation
//
      if( (tdir == trkOut ? root( f1, f2, f).success():
	   root( f2, f1, f).success()) ) {
	if(f >= fmin && f <= fmax){
//
//  Update the surface point
//
	  p = traj.position( f );  
	  d2 = surf.normalTo( p, n, uv );
	  return TrkErrCode(TrkErrCode::succeed);
	} else
	  return TrkErrCode(TrkErrCode::fail,3);
      } else
	return TrkErrCode(TrkErrCode::fail,1);
    }
//
// Try again: increment track
//
    d1 = d2;
    f1 = f2;
    iter++;
  }
// assume failure to intersect means out-of-range
  return TrkErrCode(TrkErrCode::fail,2); 
}

void
Intersection::funcd(double f, double& d, double& dd) const
{
  HepPoint   p;
  Hep3Vector n, v;
  traj.getInfo( f, p, v );
  d  =  surf.normTo( p, n );
  dd =  v*n;
}

TrkErrCode
Intersection::root(const double& x1, const double& x2, double& f) const
{
// Implement safe Newton-Raphson appropriate for traj-surface intersection 
// problem. Starts at x1, so should find solution nearest x1, and also 
// find the solution immediatedly for linear case.  If N-R step is too big, 
// (due to curved surface or traj), step by the distance to the surface instead.
// set to midpoint if starting value not in range
  double  maxDis = (x2 > x1) ? x2 : x1;
  double  minDis = (x1 < x2) ? x1 : x2;
  if(f<minDis || f>maxDis)f=(x1+x2)/2.0;
  double d, dd;
  double olddx(x2-x1);
  for( unsigned j=0; j< maxit; j++)
  {
    funcd(f, d, dd);
    if( dd==0. ) {
      f += 0.1*(x2-x1);
      continue;
    }
    double dx = d/dd;
    if(fabs(dx)<fabs(olddx))
      f += dx;
    else
      // damp the oscillations
      f += dx/2.0; 
    olddx = dx;
    if( fabs(d) < minstep ) return TrkErrCode(TrkErrCode::succeed);
  }
  return TrkErrCode(TrkErrCode::fail,2);
}

DetSurface::intertype
Intersection::findProximity(double& flen,double fmin,double fmax,
			    trkDirection tdir) const {
// first, try a linear extrapolation
// get information about the trajectories position and direction at the
// start point
  double range = fmax-fmin;
  HepPoint   pos;
  Hep3Vector dir;
  traj.getInfo( flen, pos, dir );
  double theCurvature = traj.curvature(flen);
// first test to see if an intersection if possible
  DetSurface::intermode mode = tdir==trkOut ?
    DetSurface::forward : DetSurface::backward;
  double delta;
  DetSurface::intertype inter = surf.distTo(pos,dir,delta,mode);
  double angleChange = delta*theCurvature;
  if((inter == DetSurface::intersect || inter == DetSurface::localmin) &&
     delta <= range && fabs(angleChange) < Constants::pi ) {
// success!
    flen += delta;
  } else {
// failure might be due to trajectory curvature; see if that is
// appreciable compared to the range and the current orientation
// WRT the surface
    Hep3Vector norm;
    double normdist = surf.normTo(pos,norm);
// the following is the angle change needed to head back towards the surface
    double dot = norm*dir;
    double rangle = fabs(dot)< 1.0 ? fabs(asin(dot)) : Constants::pi/2.0;
    double curve = traj.curvature(flen);
    if(curve*range > rangle){
// try stepping around the trajectory
      double rstep = rangle/curve;
      double step = std::min(std::max(std::max(rstep,fabs(normdist)),range/1000.0),range);
      if(tdir==trkIn)step *= -1.0;
// this loop allows the last step to be a little beyond the range
      while(tdir==trkOut ? flen<fmax : flen>fmin){
	flen += step;
	pos = traj.position(flen);
// check if we've crossed the surface
	double newdist = surf.normTo(pos,norm);
	if( newdist*normdist<0.0){
// success! make a linear correction back in the opposite direction
// as we came
	  flen -= step*(newdist/(newdist-normdist));
	  traj.getInfo( flen, pos, dir );
	  double delta(0.0);
	  inter = (tdir==trkOut) ?
	    surf.distTo(pos,dir,delta,DetSurface::backward) :
	    surf.distTo(pos,dir,delta,DetSurface::forward);
	  flen += delta;
	  inter = DetSurface::intersect;
	  break;
	}
	normdist=newdist;
      }
    }
  }
  return inter;
}
