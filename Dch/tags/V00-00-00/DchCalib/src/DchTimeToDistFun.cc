#include "BaBar/BaBar.hh"

#include "DchCalib/DchTimeToDistFun.hh"
#include "DchCalib/DchCalibFun.hh"
#include "DchCalib/DchGenericInverse.hh"
#include "ErrLogger/ErrLog.hh"

double DchTimeToDistFun::maxDriftTime=750;        // numerical inverse will be created covering
double DchTimeToDistFun::minDriftTime=0;          // the range [minDriftTime,maxDriftTime]
double DchTimeToDistFun::maxTimeTolerance=0.5;    // with precision better than MaxTimeTolerance 

DchTimeToDistFun::DchTimeToDistFun( std::auto_ptr<DchCalibFun>& t2d,
                                    DchCalibFun* base)
        :_fun(t2d),_inverseFun(0),_base(base)
{
}

void
DchTimeToDistFun::createInverse() const
{
        assert(_inverseFun.get()==0);
        _inverseFun.reset(_fun->inverseFunction()); // first try exact inverse
        if (_inverseFun.get()==0) {
                ErrMsg(warning) << "Generating numerical inverse..." << endmsg;
                _inverseFun.reset(new DchGenericInverse(*_fun,
                                             DchTimeToDistFun::minDriftTime,
                                             DchTimeToDistFun::maxDriftTime,
                                             DchTimeToDistFun::maxTimeTolerance)); // otherwise, approximate one
        }
        assert(_inverseFun.get()!=0);
}

void
DchTimeToDistFun::setGenericInverseParameters(double tol, double tmin,double tmax)
{
        maxTimeTolerance=tol*1e9;
        minDriftTime=tmin*1e9;
        maxDriftTime=tmax*1e9;
}
