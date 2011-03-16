//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchCalibFunVisitor.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//      Class DchCalibFunVisitor:
//         nasty little class to fake dynamic_cast for DchCalibFun's
//         because Solaris doesn't have dynamic_cast...
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gerhard Raven           11/13/2000
//
// Copyright Information:
//      Copyright (C) 2000      University of California, San Diego
//
//------------------------------------------------------------------------
#ifndef DCHCALIBFUNVISITOR_HH
#define DCHCALIBFUNVISITOR_HH

class DchCalibFun;
class DchVelocityFun;

class DchCalibFunVisitor
{
public:
        virtual ~DchCalibFunVisitor();
        virtual bool accept(const DchCalibFun&) =0; // FIXME: due to a #*($)#*$ in 
                                                    //        Solaris compiler, this one
                                                    //        _must_ be implemented in all
                                                    //        derived classes -- and everyone
                                                    //        can implement it by just 
                                                    //        returning false...
        virtual bool accept(const DchVelocityFun&);
};
#endif
