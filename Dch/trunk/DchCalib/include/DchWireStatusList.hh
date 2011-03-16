//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchWireStatusList.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//      Class DchWireStatusList
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gerhard Raven           11/01/99
//
// Copyright Information:
//      Copyright (C) 1999      University of California, San Diego
//
//------------------------------------------------------------------------
#ifndef DCHWIRESTATUSLIST_HH
#define DCHWIRESTATUSLIST_HH

#include <memory>
#include "boost/shared_ptr.hpp"
#include "DchGeom/DchWireArray.hh"
class DchWireStatus;

class DchWireStatusList
{
public:
        typedef boost::shared_ptr<DchWireStatus> DchWireStatusRCPtr;
        DchWireStatusList(std::auto_ptr<DchWireArray<DchWireStatusRCPtr> >& status);
        ~DchWireStatusList();
        const DchWireStatus* getWireStatus(unsigned layer,unsigned wire) const
        {
                return  _status->isValid(layer,wire)?_status->at(layer,wire).get():0;
        }
        template <typename I>
        const DchWireStatus* getWireStatus(const I& i) const
        {   return getWireStatus(i.layer(),i.wire()); }
private:
        DchWireStatusList& operator=(const DchWireStatusList &); // NOT IMPLEMENTED
        DchWireStatusList(const DchWireStatusList &);            // NOT IMPLEMENTED
        std::auto_ptr<DchWireArray<DchWireStatusRCPtr> >  _status;

        friend bool testCdb(const DchWireStatusList*, const DchWireStatusList*);

};
#endif
