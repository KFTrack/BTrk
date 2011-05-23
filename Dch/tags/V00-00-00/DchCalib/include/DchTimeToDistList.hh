//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchTimeToDistList.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//      Class DchTimeToDistList:
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
#ifndef DCHTIMETODISTLIST_HH
#define DCHTIMETODISTLIST_HH

#include <memory>
#include "boost/shared_ptr.hpp"
#include "DchGeom/DchWireArray.hh"

class DchTimeToDist;
class DchCalibFun;
class DchTimeToDistFun;

class DchTimeToDistList 
{
public:
        typedef boost::shared_ptr<DchTimeToDist> DchTimeToDistRCPtr;
        DchTimeToDistList(std::auto_ptr<DchWireArray<DchTimeToDistRCPtr> >& t2d);
        ~DchTimeToDistList();
        const DchTimeToDist* getTimeToDist(unsigned layer,unsigned wire) const
        {   return _t2d->isValid(layer,wire)?_t2d->at(layer,wire).get():0;
        }
        template <typename I> const DchTimeToDist* getTimeToDist(const I& i) const
        {   return getTimeToDist(i.layer(),i.wire()); }

private:
        DchTimeToDistList& operator=(const DchTimeToDistList &); // NOT IMPLEMENTED
        DchTimeToDistList(const DchTimeToDistList &);            // NOT IMPLEMENTED

        std::auto_ptr<DchWireArray<DchTimeToDistRCPtr> > _t2d;

        friend bool testCdb(const DchTimeToDistList*, const DchTimeToDistList*);  

};
#endif
