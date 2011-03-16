//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchCalibFun.cc 88 2010-01-14 12:32:57Z stroili $
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
#include "BaBar/BaBar.hh"
#include "DchCalib/DchCalibFun.hh"

DchCalibFun::DchCalibFun()
{ ; }

DchCalibFun::~DchCalibFun()
{ ; }

DchCalibFun *
DchCalibFun::inverseFunction() const 
{
        return 0;
}

bool
DchCalibFun::accept(DchCalibFunVisitor&) const 
{
        return false;
}
