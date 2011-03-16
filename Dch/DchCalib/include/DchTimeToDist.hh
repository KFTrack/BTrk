//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchTimeToDist.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchTimeToDist.
//
//	This implements a 'handle' pattern for the time to distance
//	calibration: the actual computation is delegated down to the
//	different 'strategies' for the time to distance and resolution.
//	The T0's are handled inside this class, so you could also call
//	this a 'template' pattern...
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

#ifndef DCHTIMETODIST_HH
#define DCHTIMETODIST_HH

class DchCalibFun;
class DchTimeToDistFun;
#include "boost/shared_ptr.hpp"

class DchTimeToDist
{
public:
        typedef boost::shared_ptr<DchTimeToDistFun> DchTimeToDistFunRCPtr;
        typedef boost::shared_ptr<DchCalibFun> DchCalibFunRCPtr;
        DchTimeToDist(const DchTimeToDistFunRCPtr& tFun,
                      const DchCalibFunRCPtr& rFun,
                      double t0);
        DchTimeToDist(const DchTimeToDist &rhs);
        ~DchTimeToDist();

        double timeToDist(double time,
                          int ambiguity=0,        // left or right?
                          double entranceAngle=0,
                          double dipAngle=0,
                          double z=0,
                          double Q=0) const;

        double resolution(double distance,
                          int ambiguity=0,
                          double entranceAngle=0,
                          double dipAngle=0,
                          double z=0,
                          double Q=0) const;

        double distToTime(double distance,
                          int ambiguity,
                          double entranceAngle,
                          double dipAngle,
                          double z,
                          double Q=0 ) const;

protected:
// these last two are specific to calibration code. Don't use them unless
// you have a _very_ good reason to!
        friend class DchTimeToDistFunAccessor;
        const DchTimeToDistFun* timeToDistFun() const { return _fun.get();}
        const DchCalibFun* resolutionFun() const { return _rFun.get();}
        double t0() const { return _t0; }

private:
        DchTimeToDistFunRCPtr _fun;
        DchCalibFunRCPtr      _rFun;  // resolution function
        double                _t0;

        DchTimeToDist();            // not implemented (on purpose!!!)
        DchTimeToDist& operator =(const DchTimeToDist&);

        // speed-up hack
        unsigned _tPar,_rPar;

        friend bool testCdb(const DchTimeToDist*, const DchTimeToDist*);
};

#endif // DCHTIMETODIST_HH
