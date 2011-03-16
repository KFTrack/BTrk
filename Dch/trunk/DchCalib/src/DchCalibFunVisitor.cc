//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchCalibFunVisitor.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//      Class DchCalibFunVisitor:
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
#include "BaBar/BaBar.hh"
#include "DchCalib/DchCalibFunVisitor.hh"

DchCalibFunVisitor::~DchCalibFunVisitor()
{ ; }

bool
DchCalibFunVisitor::accept(const DchCalibFun&)
{
        return false;
}

bool
DchCalibFunVisitor::accept(const DchVelocityFun&)
{
        return false;
}
