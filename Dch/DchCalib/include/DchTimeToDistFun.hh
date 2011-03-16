#ifndef DCHTIMETODISTFUN
#define DCHTIMETODISTFUN

#include "DchCalib/DchCalibFun.hh"
#include <memory>

class DchTimeToDistFun
{
public:
        DchTimeToDistFun(std::auto_ptr<DchCalibFun>& t2d,
                         DchCalibFun* base);

        double value(const HepVector& x) const { return _fun->value(x); }
        unsigned numberOfArguments() const     { return _fun->numberOfArguments(); }
        const char *channelName() const        { return _fun->channelName(); }
        double inverseValue(const HepVector& x) const { if (_inverseFun.get()==0) createInverse();
                                                        return _inverseFun->value(x); }
        const DchCalibFun *base() const        { return _base; }
        const DchCalibFun *fun() const         { return _fun.get(); }

        static void setGenericInverseParameters(double tol,double tmin, double tmax);
                                                // when creating a numerical inverse, make
                                                // sure it covers the range [tmin,tmax], and
                                                // the max deviation is no larger than 'tol'
                                                // NOTE: units are SECONDS
private:
        static double maxTimeTolerance;
        static double maxDriftTime,minDriftTime;

        std::auto_ptr<DchCalibFun> _fun;
        mutable std::auto_ptr<DchCalibFun> _inverseFun;
        DchCalibFun *_base;

        void createInverse() const;

        friend bool testCdb(const DchTimeToDistFun*, const DchTimeToDistFun*); 

};

#endif
