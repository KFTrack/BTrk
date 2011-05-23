//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchGenericInverse.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//      Class DchGenericInverse:
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
#ifndef DCHGENERICINVERSE_HH
#define DCHGENERICINVERSE_HH

#include <assert.h>
#include "DchCalib/DchCalibFun.hh"
class DchGenericInverseLUT;


class DchGenericInverse: public DchCalibFun 
{
public:
        DchGenericInverse(const DchCalibFun& ,double mintime, double maxtime, 
                          double tolerance);
        ~DchGenericInverse();

        double value(const HepVector &) const;
        HepMatrix gradientWrtParameters(const HepVector& arguments)  const;
        DchCalibFun* inverseFunction() const;
        DchCalibFun* clone() const;
        unsigned numberOfArguments() const;
        unsigned numberOfParameters() const;

        const char *fieldName(unsigned i) const;
        const char *fieldUnit(unsigned i) const;
        double      fieldValue(unsigned i) const;
        const char* channelName() const ;

        DchGenericInverse(const DchGenericInverse &);
private:
        DchGenericInverse& operator=(const DchGenericInverse &); // NOT IMPLEMENTED

        const DchCalibFun& _fun;
        DchGenericInverseLUT *_tbl,*_tbll,*_tblr;
        double _tolerance;
        unsigned _nArg;

        const DchGenericInverseLUT* amb2tbl(double i) const
        { return i==0?_tbl:i<0?_tbll:_tblr; }
        DchGenericInverseLUT*& amb2tbl(double i) 
        { return i==0?_tbl:i<0?_tbll:_tblr; }
};
#endif
