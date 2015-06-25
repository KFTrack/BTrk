//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHotList.cc,v 1.22 2004/08/06 06:31:42 bartoldu Exp $
//
// Description:
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkHotList.hh"
#include "BTrk/TrkBase/TrkHitOnTrk.hh"
#include "BTrk/TrkBase/TrkView.hh"
#include "BTrk/BaBar/BbrCollectionUtils.hh"
#include <iostream>
using std::endl;
using std::ostream;

TrkHotList::TrkHotList()
{
}

TrkHotList::~TrkHotList()
{
}

void
TrkHotList::print(ostream& o) const
{
       o << " hitCapable: " << (hitCapable()?"yes":"no")
         << " nActive: " << nActive()
         << " nHit: " << nHit()
         << " startFoundRange: " <<startFoundRange()
         << " endFoundRange: " << endFoundRange();
}

void
TrkHotList::printAll(ostream &o) const
{
        print(o); o << "\n";
        TrkHotList::hot_iterator i= begin();
        while (i!=end()) {
                i++->print(o); o << endl;
        }
}


TrkHotList*
TrkHotList::resetParent(TrkBase::Functors::setParent f)
{
        std::for_each(begin(),end(),f);
        return this;
}

void
TrkHotList::sort()
{
        std::sort(hotlist().begin(),
                  hotlist().end(),
                  babar::Collection::PtrLess());
}
