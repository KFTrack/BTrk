//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchGenericInverse.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//      Class DchGenericInverse
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gerhard Raven           8/05/98
//
// Copyright Information:
//      Copyright (C) 1998      University of California, San Diego
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include "DchCalib/DchGenericInverse.hh"
#include "DchCalib/DchGenericInverseLUT.hh"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "ErrLogger/ErrLog.hh"

DchGenericInverse::DchGenericInverse(const DchCalibFun& fun, 
                                     double mintime, double maxtime, 
                                     double tolerance)
        :_fun(fun),_tbl(0),_tbll(0),_tblr(0),_tolerance(tolerance),
         _nArg(fun.numberOfArguments())
{
        assert(maxtime>mintime);

        ErrMsg(routine) << "generating numerical inverse with tolerance " 
                        << tolerance << endmsg;
        if (_nArg==1) {
                HepVector x(1);
                amb2tbl(0)= new DchGenericInverseLUT(fun,x,mintime,maxtime,
                                                     tolerance);
                ErrMsg(routine) << "generated single LUT with " 
                                << amb2tbl(0)->entries() 
                                << " entries" << endmsg;
        } else if (_nArg==2) {
                HepVector x(2);
                static const double amb[3]={-1,0,1};
                for (unsigned i=0;i<3;++i) {
                        x[1]=amb[i];
                        amb2tbl(amb[i]) = 
                                new DchGenericInverseLUT(fun,x,mintime,maxtime,
                                                         tolerance);
                }
                ErrMsg(routine) << "generated three LUTs with " 
                                << amb2tbl(0)->entries() << "," 
                                << amb2tbl(-1)->entries() << "," 
                                << amb2tbl(1)->entries() 
                                << " entries" << endmsg;
        } else if (_nArg==3) {
                HepVector x(3);
		x[2]=0.0; //calculate inverse assuming entrance angle = 0.0
                static const double amb[3]={-1,0,1};
                for (unsigned i=0;i<3;++i) {
                        x[1]=amb[i];
                        amb2tbl(amb[i])= 
                                new DchGenericInverseLUT(fun,x,mintime,maxtime,
                                                         tolerance);
                }
                ErrMsg(routine) << "generated three LUTs with " 
                                << amb2tbl(0)->entries() << "," 
                                << amb2tbl(-1)->entries() << "," 
                                << amb2tbl(1)->entries() 
                                << " entries" << endmsg;
        } else {
                ErrMsg(fatal) << " Generic inverse not yet configured for functions"
                              << " with more than 3 variables... " << endmsg;
        }
}

DchGenericInverse::DchGenericInverse(const DchGenericInverse &other)
        : _fun(other._fun),
          _tbl(new DchGenericInverseLUT(*other._tbl)),
          _tbll(new DchGenericInverseLUT(*other._tbll)),
          _tblr(new DchGenericInverseLUT(*other._tblr)),
          _tolerance(other._tolerance),
          _nArg(other._nArg)
{
}

DchGenericInverse::~DchGenericInverse()
{ 
        delete _tbl;
        delete _tbll; 
        delete _tblr;
}

HepMatrix 
DchGenericInverse::gradientWrtParameters(const HepVector& arguments)  const
{
        ErrMsg(fatal) << " You're asking for gradients of a generic inverse;"
                      << " This hasn't been required sofar, and therefore it"
                      << " hasn't been implemented yet..." << endmsg;
        return HepMatrix();
}

DchCalibFun*
DchGenericInverse::inverseFunction() const 
{
        return _fun.clone();  // return a clone of the original function...
}

DchCalibFun*
DchGenericInverse::clone() const
{ 
        return new DchGenericInverse(*this);
}

unsigned 
DchGenericInverse::numberOfArguments() const
{
        return _fun.numberOfArguments();
}

unsigned 
DchGenericInverse::numberOfParameters() const 
{
        unsigned npar=0;
        npar=amb2tbl(0)->entries();
        if (_nArg>1) {
                npar+=amb2tbl(-1)->entries();
                npar+=amb2tbl(+1)->entries();
        }
        return 2*npar;
}

const char *
DchGenericInverse::fieldName(unsigned i) const 
{
        return i%2?"t(i)":"d(i)";
}

const char *
DchGenericInverse::fieldUnit(unsigned i) const
{
        return i%2?"ns":"cm";
}

double      
DchGenericInverse::fieldValue(unsigned i) const
{
        unsigned n=2*amb2tbl(0)->entries();
        if (i<n) return amb2tbl(0)->fieldValue(i);
        if (_nArg>1) {
                i-=n; 
                n=2*amb2tbl(-1)->entries();
                if (i<n) return amb2tbl(-1)->fieldValue(i);
                i-=n;
                n=2*amb2tbl(+1)->entries();
                if (i<n) return amb2tbl(+1)->fieldValue(i);
        }
        ErrMsg(error) << " not that many fields..." << endmsg;
        return 0;
}

const char * 
DchGenericInverse::channelName() const 
{
        return "DchGenericInverse";
}

double 
DchGenericInverse::value(const HepVector &arg) const
{
        assert(arg.num_row()==_nArg);
        static HepVector t1(1),t2(2),t3(3); HepVector *pt=0;
        double amb=0;
        if (_nArg==1) {
                pt=&t1;
        } else if (_nArg==2) {
                pt=&t2; amb=arg[1];
                t2[1]=arg[1]; 
        } else if (_nArg==3) {
                pt=&t3; amb=arg[1];
                t3[1]=arg[1]; 
                t3[2]=arg[2]; 
        }
        HepVector &t=*pt;
        const DchGenericInverseLUT &lut=*amb2tbl(amb);

        double d=arg[0];
        t[0]=lut.value(d);
        // check quality...
        if ( fabs(_fun.value(t)-d)>_tolerance ) {
                ErrMsg(warning) 
                        << "numerical inverse exceeds tolerance: fabs("
                        << _fun.value(t) << " - " << d << ") > " 
                        << _tolerance << endmsg;
        } 
        return t[0];
}
