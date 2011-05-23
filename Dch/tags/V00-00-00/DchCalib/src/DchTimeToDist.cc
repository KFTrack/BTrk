//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchTimeToDist.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchTimeToDist
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	G. Raven
//	
//
// Copyright Information:
//	Copyright (C) 1998	University of California, San Diego
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------
#include "DchCalib/DchTimeToDist.hh"
#include "DchCalib/DchTimeToDistFun.hh"
#include "DchCalib/DchCalibFun.hh"
#include "DchCalib/DchGenericInverse.hh"
#include "DchGeom/DchDetector.hh"
#include <assert.h>
#include <math.h>
#include "CLHEP/Matrix/Vector.h"
#include "ErrLogger/ErrLog.hh"

static const char rscid[] = "$Id: DchTimeToDist.cc 88 2010-01-14 12:32:57Z stroili $";
static const double wireRadius=0.0010;
static const double maxdist=1.5;

DchTimeToDist::DchTimeToDist(const DchTimeToDistFunRCPtr& fun,
                             const DchCalibFunRCPtr& rFun,
                             double t0)
       :
       _fun(fun),_rFun(rFun), _t0(t0),
       _tPar(fun->numberOfArguments()),
       _rPar(rFun->numberOfArguments())
{
        assert(_tPar>0&&_tPar<6);
        assert(_rPar>0&&_rPar<6);
}

DchTimeToDist::DchTimeToDist(const DchTimeToDist &rhs)
        :_fun(rhs._fun),_rFun(rhs._rFun),_t0(rhs._t0),
         _tPar(rhs._tPar),_rPar(rhs._rPar)
{
}


DchTimeToDist::~DchTimeToDist()
{
}

double
DchTimeToDist::timeToDist(double time,
                          int    ambiguity,
                          double entranceAngle,
                          double dipAngle,
                          double z,
                          double Q) const
{
        // Use static vectors to avoid overhead of HepVector
        //   creation and destruction.  Ugly but effective.--sfs
        static int dmutex=0;
        static HepVector p1(1),p2(2),p3(3),p4(4),p5(5),p6(6);
        static HepVector* pptr[]={&p1,&p2,&p3,&p4,&p5,&p6};
        HepVector& p = *pptr[_tPar-1];


        if (++dmutex!=1) ErrMsg(fatal) << "re-entering non-reentrant code" << endmsg;
        switch (_tPar) {  // simple version of Duff's device ;-)
                case 6: p[5]=Q;
                case 5: p[4]=z;
                case 4: p[3]=dipAngle;
                case 3: p[2]=entranceAngle;
                case 2: p[1]=ambiguity;
                case 1: p[0]=time-_t0;
        }
        double dist=_fun->value(p);
        --dmutex;

        return dist<wireRadius?wireRadius:(dist>maxdist?maxdist:dist);
}

double
DchTimeToDist::resolution(double distance,
                          int    ambiguity,
                          double entranceAngle,
                          double dipAngle, double z, double Q) const
{
        // Use static vectors to avoid overhead of HepVector
        //   creation and destruction.  Ugly but effective.--sfs
        static int rmutex=0;
        static HepVector p1(1),p2(2),p3(3),p4(4),p5(5),p6(6);
        static HepVector* pptr[]={&p1,&p2,&p3,&p4,&p5,&p6};
        HepVector& p = *pptr[_rPar-1];

        if (++rmutex!=1) ErrMsg(fatal) << "re-entering non-reentrant code" << endmsg;
        switch (_rPar) {  // simple version of Duff's device ;-)
                case 6: p[5]=Q;
                case 5: p[4]=z;
                case 4: p[3]=dipAngle;
                case 3: p[2]=entranceAngle;
                case 2: p[1]=ambiguity;
                case 1: p[0]=fabs(distance);
        }
        double r=_rFun->value(p);
        --rmutex;
        if (r<=0.0030) {
                ErrMsg(warning) << "computed resolution = " << r << endmsg;
        }
        return r;
}

double
DchTimeToDist::distToTime(double distance, int ambiguity,
                          double entranceAngle,
                          double dipAngle, double z, double Q) const
{
        // Use static vectors to avoid overhead of HepVector
        //   creation and destruction.  Ugly but effective.--sfs
        static int tmutex=0;
        static HepVector p1(1),p2(2),p3(3),p4(4),p5(5),p6(6);
        static HepVector* pptr[6]={&p1,&p2,&p3,&p4,&p5,&p6};

        HepVector& p = *pptr[_tPar-1];

        if (++tmutex!=1) ErrMsg(fatal) << "re-entering non-reentrant code" << endmsg;
        switch (_tPar) {
                case 6: p[5]=Q;
                case 5: p[4]=z;
                case 4: p[3]=dipAngle;
                case 3: p[2]=entranceAngle;
                case 2: p[1]=ambiguity;
                case 1: p[0]=fabs(distance);
        }
        double t=_fun->inverseValue(p)+_t0;
        --tmutex;
        return t;
}
