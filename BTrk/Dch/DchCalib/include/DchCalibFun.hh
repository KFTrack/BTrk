//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchCalibFun.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//      Class DchCalibFun:
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
#ifndef DCHCALIBFUN_HH
#define DCHCALIBFUN_HH

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
using namespace CLHEP;

#include <iosfwd>
class DchCalibFunVisitor;

class DchCalibFun
{
public:
        DchCalibFun();
        virtual ~DchCalibFun();

        virtual double value(const HepVector& arguments) const=0;
        virtual HepMatrix gradientWrtParameters(const HepVector& arguments)  const=0;
        virtual DchCalibFun *inverseFunction() const; // if your function is invertible, 
                                                      // this is where you return the inverse function!
        virtual DchCalibFun *clone() const =0;
        virtual unsigned numberOfArguments() const =0;
        virtual unsigned numberOfParameters() const =0;

        virtual bool accept(DchCalibFunVisitor& ) const;

        virtual const char *fieldName(unsigned i) const=0;
        virtual const char *fieldUnit(unsigned i) const=0;
        virtual double      fieldValue(unsigned i) const=0;
        virtual const char* channelName() const =0;
};
#endif
