//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchDedxCalibList.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//      Class DchDedxCalibList:
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Fergus Wilson           02-FEB-1999
//
// Changes:
//      Jean Roy     7-JAN-2000
//          Replace the multiple calibration functions by a single
//          combined function.
//
// Copyright Information:
//      Copyright (C) 1999      University of California, San Diego
//      Copyright (C) 2000      University of Colorado, Boulder
//
//------------------------------------------------------------------------
#ifndef DCHDEDXCALIBLIST_HH
#define DCHDEDXCALIBLIST_HH

#include <assert.h>
#include <memory>
#include <vector>

#include <iostream>
#include "boost/shared_ptr.hpp"
#include "DchGeom/DchWireArray.hh"
#include "DchCalib/DchDedxCalib.hh"

class DchDedxCalib;
class DchCalibFun;

/**
  * Class for managing all the dE/dx calibration objects
  * @author Fergus Wilson
  * @version $Id: DchDedxCalibList.hh 88 2010-01-14 12:32:57Z stroili $
  **/
class DchDedxCalibList 
{
public:
     typedef boost::shared_ptr<DchCalibFun> DchCalibFunRCPtr;
     typedef boost::shared_ptr<DchDedxCalib> DchDedxCalibRCPtr;
  /**
    * @param gainFactors   Table of pointers to gainFactor correction 
    * numbers accessed by (layer,wire). Although the gain is single number, 
    * it has the same shape as the following tables.
    * @param doca          List of pointers to doca correction objects
    * @param docamap       Table of pointers to doca correction objects accesed
    *                      by (layer,wire)
    * @param entranceAngle List of pointers to Entrance Angle correction objects
    * @param entranceAnglemap Table of pointers to Entrance Angle correction objects accesed
    *                      by (layer,wire)
    * @param dipAngle      List of pointers to Dip Angle correction objects
    * @param dipAnglemap   Table of pointers to Dip Angle correction objects accesed
    *                      by (layer,wire)
    * @param z             List of pointers to z correction objects
    * @param zmap          Table of pointers to z correction objects accesed
    *                      by (layer,wire)
    */
  DchDedxCalibList(double globalGasGain,
                   std::auto_ptr<DchWireArray<double> >& gainFactors,
                   std::auto_ptr<DchWireArray<DchCalibFunRCPtr> >& corrFunMap);

  /**
    * Destructor
    */
  ~DchDedxCalibList();

  /**
    * function to check if object already created and, if not, create a
    * new DchDedxCalib object and add it to the DchDedxCalibList; finally
    * return pointer to DchDedxCalib.
    * @return returns pointer to DchDedxCalib object for this (layer,wire)
    * @param layer Layer number [1,40]
    * @param wire  Wire number [1,256]
    */
  const DchDedxCalib* getDedxCalib(unsigned layer,unsigned wire) const;

   template <typename I>
   const DchDedxCalib* getDedxCalib(const I& i) const
   {   return getDedxCalib(i.layer(),i.wire()); }

private: 
  DchDedxCalibList& operator=(const DchDedxCalibList&); // NOT IMPLEMENTED
  DchDedxCalibList(const DchDedxCalibList&);            // NOT IMPLEMENTED

  double _globalGasGain;
  // A 2-D table of gain factors accessed by (layer,wire).
  std::auto_ptr<DchWireArray<double> > _gainFactors;

  // A list of pointers to calibration objects for the overall correction,
  // organized by (layer,wire).
  // This function is the combination of several functions and is produced
  // by the class DchDedxCorrFun.
  // This list owns the objects so that when this list is deleted the
  // objects are also deleted. 
  std::auto_ptr<DchWireArray<DchCalibFunRCPtr> > _corrFun;

  // Table of calibration objects. It is just a speedup aid to check
  // if an object has already been created and returns it to the caller.
  mutable DchWireArray<DchDedxCalibRCPtr>  _cache;

  friend bool testCdb(const DchDedxCalibList*, const DchDedxCalibList*);

};
#endif
