//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPocaXY.cc,v 1.11 2006/03/25 15:15:56 brownd Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner, largely taken from Art Snyder
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "ErrLogger/ErrLog.hh"
#include "TrkBase/TrkPocaXY.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/HelixParams.hh"
#include "TrkBase/HelixTraj.hh"
#include "TrajGeom/TrkLineTraj.hh"
#include "CLHEP/Vector/ThreeVector.h"
using std::endl;

TrkPocaXY::~TrkPocaXY()
{}

double
TrkPocaXY::doca() const {
  return _docaxy;
}

TrkPocaXY*
TrkPocaXY::clone() const {
  return new TrkPocaXY(*this);
}


TrkPocaXY::TrkPocaXY(const TrkPocaXY& other) :
  TrkPocaBase(other),_docaxy(other._docaxy)
{}

TrkPocaXY::TrkPocaXY(const Trajectory& traj, double flt, 
                 const HepPoint& pt, double prec) 
  : TrkPocaBase(flt,prec) , _docaxy(-9999.0) {
// construct a line trajectory through the given point
  Hep3Vector zaxis(0.0,0.0,1.0);
// length is set to 200cm (shouldn't matter)
  TrkLineTraj zline(pt, zaxis, 200.0);
// create a 2-traj poca with this
  double zlen(pt.z());
  TrkPoca zlinepoca(traj,flt,zline,zlen,prec);
// transfer over the data
  _status = zlinepoca.status();
  if(status().success()){
    _flt1 = zlinepoca.flt1();
    _flt2 = zlinepoca.flt2();
    _docaxy  = zlinepoca.doca(); // doca should be perpendicular to the zline so 2d by construction
  }
}

TrkPocaXY::TrkPocaXY(const Trajectory& traj1,double flt1,
		     const Trajectory& traj2,double flt2,
		     double prec) 
  : TrkPocaBase(flt1,flt2,prec) , _docaxy(-9999.0) {

  // this traj-traj version starts by approximating the trejectories as
  // either a line or a circle and treats separately the line-line, 
  // circle-circle, and line-circle cases.

  _flt1 = flt1;
  _flt2 = flt2;
  double delta1(10*prec);
  double delta2(10*prec);
  unsigned niter(0);
  static const unsigned maxiter(20);
  while(niter < maxiter && 
	( delta1 > prec || delta2 > prec ) ){
    // get positions and directions
    HepPoint  pos1 = traj1.position(_flt1);
    HepPoint pos1b;
    if(niter == 0)  pos1b = pos1;
    HepPoint  pos2 = traj2.position(_flt2);
    Hep3Vector dir1 = traj1.direction(_flt1);
    Hep3Vector dir2 = traj2.direction(_flt2); 
    Hep3Vector dd1 = traj1.delDirect(_flt1);
    Hep3Vector dd2 = traj2.delDirect(_flt2);
    double curv1 = traj1.curvature(_flt1); 
    double curv2 = traj2.curvature(_flt2); 
    double m1, m2, q1, q2;
    double r1, r2, xc1,xc2,yc1,yc2,sr1,sr2;
    if(curv1 == 0) {
      // convert to m*x+q
      if(dir1[0] == 0){
        m1 = 1.0e12;
      }else{
        m1 = dir1[1]/dir1[0];
      }
      q1 = pos1.y()-pos1.x()*m1;
    }else{
    // get circle parameters
      r1 = (1- dir1[2]*dir1[2])/curv1;
      sr1=r1;
      if(dir1[0]*dd1[1] - dir1[1]*dd1[0] < 0) sr1 = -r1;
      double cosphi1 = dir1[0]/sqrt(dir1[0]*dir1[0]+dir1[1]*dir1[1]);
      double sinphi1 = cosphi1 * dir1[1]/dir1[0];
      xc1 = pos1.x() - sr1 * sinphi1;
      yc1 = pos1.y() + sr1 * cosphi1;
    }
    if(curv2 == 0) {
      if(dir2[0] == 0){
        m2 = 1.0e12;
      }else{ 
        m2 = dir2[1]/dir2[0];
      }
       q2 = pos2.y()-pos2.x()*m2;
    }else{
      r2 = (1-dir2[2]*dir2[2])/curv2;
      sr2 = r2;
      if(dir2[0]*dd2[1] - dir2[1]*dd2[0] < 0) sr2 = -r2;
      double cosphi2 = dir2[0]/sqrt(dir2[0]*dir2[0]+dir2[1]*dir2[1]);
      double sinphi2 = cosphi2 * dir2[1]/dir2[0];
      xc2 = pos2.x() - sr2 * sinphi2;
      yc2 = pos2.y() + sr2 * cosphi2;
    }
    double xint, yint, xint1, xint2, yint1, yint2, absdoca;
    _docaxy = 0;
    _status.setSuccess(3);

    // First the line-line case

    if(curv1==0 && curv2==0){
    // do the intersection in 2d
      interTwoLines(m1, q1, m2, q2, xint, yint);
    }

    //  next do the two circle case

    if(curv1 !=0 && curv2 !=0){
      //There are four cases
      double cdist = sqrt((xc1-xc2)*(xc1-xc2)+(yc1-yc2)*(yc1-yc2));
      // a - coincident centers
      if (fabs(cdist) < 1.e-12 ) {      
       //the algorithm fails because the points have all the same distance
        _status.setFailure(12, 
                  "TrkPocaXY:: the two circles have the same center...");
         return;
      }
     // b - Actual intersection
      if ( (fabs(r1-r2) < cdist) && (cdist < r1+r2 ) ) {      
        interTwoCircles(xc1,yc1,r1,xc2,yc2,r2,xint1,yint1,xint2,yint2);
        double dist1 = (xint1-pos1b.x())*(xint1-pos1b.x()) + 
                       (yint1-pos1b.y())*(yint1-pos1b.y());
        double dist2 = (xint2-pos1b.x())*(xint2-pos1b.x()) + 
                       (yint2-pos1b.y())*(yint2-pos1b.y());

        //choose closest to begining of track1
        if(dist1<dist2){
	  xint = xint1;
          yint = yint1;
        } else {
          xint = xint2; 
          yint = yint2;
        }    
      }

    // c - nested circles and d - separated circles

      if ( (fabs(r1-r2) > cdist) ||  // nested circles
         (    cdist > (r1+r2) )) {   // separated circles
        // line going through the centers of the two circles

        double m = (yc1-yc2)/(xc1-xc2);  // y = m*x+q
        double q = yc1 - m*xc1;

        // intersection points between the line and the two circles

        double xint1, yint1, xint2, yint2, zOfDCrossT1;

        interLineCircle(m, q, xc1, yc1, r1, xint1, yint1, xint2, yint2);

        double xint3, yint3, xint4, yint4 ;
        interLineCircle(m, q, xc2, yc2, r2, xint3, yint3, xint4, yint4);
        if (fabs(r1-r2) > cdist) { // nested circles
          absdoca = fabs(r1-r2)-cdist;
          ErrMsg(debugging) << " doing nested circles in TrkPocaXY " << endmsg;

          double dist1_3 = pow((xint1-xint3),2.) + pow((yint1-yint3),2.);
          double dist2_4 = pow((xint2-xint4),2.) + pow((yint2-yint4),2.);

          if(dist1_3<dist2_4) {
            xint = 0.5*(xint1+xint3);
            yint = 0.5*(yint1+yint3);
            zOfDCrossT1 = (xint3-xint1)*dir1[1]-(yint3-yint1)*dir1[0];
          } else {
            xint = 0.5*(xint2+xint4);
            yint = 0.5*(yint2+yint4);
            zOfDCrossT1 = (xint4-xint2)*dir1[1]-(yint4-yint2)*dir1[0];
	  }  
          _docaxy = absdoca;
          if( zOfDCrossT1 > 0) _docaxy = -absdoca;   
	}
        if( cdist > (r1+r2) ) { //separated circles
          absdoca = cdist - (r1+r2);
	  ErrMsg(debugging)<<" doing separated circles in TrkPocaXY"<< endmsg;

          double dist2_3 = pow((xint2-xint3),2.) + pow((yint2-yint3),2.);
          double dist1_4 = pow((xint1-xint4),2.) + pow((yint1-yint4),2.);
          if (dist2_3<dist1_4){
            xint = 0.5*(xint2+xint3);
            yint = 0.5*(yint2+yint3);
            zOfDCrossT1 = (xint3-xint1)*dir1[1]-(yint3-yint1)*dir1[0];
          } else {
            xint = 0.5*(xint1+xint4);
            yint = 0.5*(yint1+yint4);
            zOfDCrossT1 = (xint4-xint2)*dir1[1]-(yint4-yint2)*dir1[0];
	  }
          _docaxy = absdoca;
          if( zOfDCrossT1 > 0) _docaxy = -absdoca;   
	}
      }
    }

    // Now do the line-circle case

    if((curv1 == 0 && curv2 !=0) || (curv1 != 0 && curv2 == 0)) {
  // distance between the line and the circle center
      HepVector dirT;
      double m,q,r,xc,yc, zOfDCrossT1;
      if(curv1==0) {m=m1; q=q1; r=r2; xc=xc2; yc=yc2; dirT=dir2;}
      else          {m=m2; q=q2; r=r1; xc=xc1; yc=yc1; dirT=dir1;}

      double dist = fabs((m*xc-yc+q)/sqrt(1+m*m));
      if(dist <= r) {

        // the intersection points

        interLineCircle(m,q,xc,yc,r, xint1, yint1, xint2, yint2);
        double dist1 = (xint1-pos1b.x())*(xint1-pos1b.x()) + 
                       (yint1-pos1b.y())*(yint1-pos1b.y());
        double dist2 = (xint2-pos1b.x())*(xint2-pos1b.x()) + 
                       (yint2-pos1b.y())*(yint2-pos1b.y());
        //choose closest to the beginning of track 1
        if(dist1<dist2){
  	  xint = xint1;
          yint = yint1;
        } else {
          xint = xint2;
          yint = yint2;
        }    
      } else { // no intersection points
	ErrMsg(debugging)<<" doing separated line-circle in TrkPocaXY"<<endmsg;

        // line going through the circle center and perpendicular to traj1

        double mperp = -1./m;
        double qperp = yc - mperp*xc;

        // intersection between this line and the two trajectories

        interLineCircle(mperp, qperp, xc, yc, r, xint1,yint1,xint2,yint2);

        double xint3,yint3;
        interTwoLines(m, q, mperp, qperp, xint3, yint3);

        double dist1_3 = pow((xint1-xint3),2.) + pow((yint1-yint3),2.);
        double dist2_3 = pow((xint2-xint3),2.) + pow((yint2-yint3),2.);
        if (dist1_3<dist2_3) {
          xint =  0.5*(xint3 + xint1);
          yint =  0.5*(yint3 + yint1);
          absdoca = sqrt((xint3-xint1)*(xint3-xint1)+
                         (yint3-yint1)*(yint3-yint1));
          zOfDCrossT1 = (xint3-xint1)*dirT[1]-(yint3-yint1)*dirT[0];
        }
        else {
          xint = 0.5*(xint3 + xint2);
          yint = 0.5*(yint3 + yint2);
          absdoca = sqrt((xint3-xint2)*(xint3-xint2)+
                         (yint3-yint2)*(yint3-yint2));
          zOfDCrossT1 = (xint3-xint2)*dirT[1]-(yint3-yint2)*dirT[0];
	}
         _docaxy = absdoca;
        if( zOfDCrossT1 > 0) _docaxy = -absdoca;   
      }
    }

    // get the flight lengths for the intersection point

    HepPoint intpt( xint, yint, 0.0);
    TrkPocaXY poca1(traj1,_flt1,intpt,prec);
    TrkPocaXY poca2(traj2,_flt2,intpt,prec);
    _status = poca2.status();
    if(poca1.status().success() && poca2.status().success()) {
       delta1 = fabs(_flt1 - poca1.flt1());
       delta2 = fabs(_flt2 - poca2.flt1());
       _flt1 = poca1.flt1();
       _flt2 = poca2.flt1(); 
    }else{
      ErrMsg(warning) << " poca fialure " << endmsg;
      if(poca1.status().failure()) _status=poca1.status();
      break;
    }
    niter++;
  }
  ErrMsg(debugging) <<"TrkPocaXY iterations = " << niter-1 
       << "       flight lengths " << _flt1 <<" " << _flt2 << endl 
       << "       deltas " << delta1 <<"  " << delta2 << endmsg;
  if(niter == maxiter){
    ErrMsg(routine) << "TrkPocaXY:: did not converge" << endmsg;
    _status.setFailure(13, "TrkPocaXY:: did not converge"); 
  }
}
//------------------------------------------------------------------------
void
TrkPocaXY::interLineCircle(const double& m, const double& q,
                    const double& xc, const double& yc, const double& r,
                    double& xint1,  double& yint1,
                    double& xint2,  double& yint2)
 //-------------------------------------------------------------------------
{

 double alpha = 1 + m*m;

 double beta = -xc +m*(q-yc);

 double gamma = xc*xc + (q-yc)*(q-yc) -r*r;

 double Delta = beta*beta - alpha*gamma;

 if (Delta < 0) {

    _status.setFailure(14,
              "TrkPocaXY:: the line and the circle heve no intersections...");
    return;

  }
  else if (fabs(Delta) <1.e-12) {

  xint1 = -beta/alpha;
  xint2 = xint1;

  }
  else {

    double xPlus = -beta/alpha + sqrt(beta*beta - alpha*gamma)/alpha;
    double xMinus = -beta/alpha - sqrt(beta*beta - alpha*gamma)/alpha;

    if (xPlus > xMinus) {
      xint1 = xMinus;
      xint2 = xPlus;
    }
    else {
      xint1 = xPlus;
      xint2 = xMinus;
    }
  }
 yint1 = m*xint1 + q;
 yint2 = m*xint2 + q;

 return;
}
 //------------------------------------------------------------------------
void
TrkPocaXY::interTwoCircles
(const double& xc1, const double& yc1, const double& r1,
 const double& xc2, const double& yc2, const double& r2,
 double& xint1,  double& yint1, double& xint2,  double& yint2)
  //-------------------------------------------------------------------------
{

  double A = (xc1*xc1 + yc1*yc1 - r1*r1) - (xc2*xc2 + yc2*yc2 - r2*r2);

  double B = -xc1 + xc2;

  double C = -yc1 + yc2;

  double alpha = 1 + (B*B)/(C*C);

  double beta = -xc1 + B/C*(yc1+A/(2*C));

  double gamma = xc1*xc1 + (yc1+A/(2*C))*(yc1+A/(2*C)) - r1*r1;

  double Delta = beta*beta - alpha*gamma;

  if (Delta < 0) {

    _status.setFailure(14, "TrkPocaXY:: the two circles have no intersections..\
.");
    return;

  }
  else if (fabs(Delta) <1.e-12) {

    xint1 = -beta/alpha;
    xint2 = xint1;

  }
  else {
    double xPlus = -beta/alpha + sqrt(beta*beta - alpha*gamma)/alpha;
    double xMinus = -beta/alpha - sqrt(beta*beta - alpha*gamma)/alpha;

    if (xPlus > xMinus) {
      xint1 = xMinus;
      xint2 = xPlus;
    }
    else {
      xint1 = xPlus;
      xint2 = xMinus;
    }

  }

  yint1 = -(A+2*B*xint1)/(2*C);
  yint2 = -(A+2*B*xint2)/(2*C);


  return;
}

//------------------------------------------------------------------------
void	
TrkPocaXY::interTwoLines
(const double& m1, const double& q1, const double& m2, const double& q2,
 double& xint,  double& yint)
 //-------------------------------------------------------------------------
{

  if (fabs(m1-m2) < 1.e-12) {  // parallel lines

    //the algorithm fails because the points have all the same distance

    _status.setFailure(13, "TrkPocaXY:: the two lines are parallel...");
    return;
  }
  else { // the lines have an intersection point

    xint = (q2-q1)/(m1-m2);
    yint = m1*xint + q1;
  }

return;
}
