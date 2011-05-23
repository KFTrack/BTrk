//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchDedxCalibList.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchDedxCalibList
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Fergus Wilson           02-FEB-1999
//
// Changes:
//      Jean Roy   7-JAN-2000
//          Replace the multiple calibration functions by a single
//          combined function.
//
//	Michael Kelsey 10-May-2004
//	    See DchTimeToDistList:  modify access function to "skip" channels
//	    with null pointer, rather than asserting.  This supports using
//	    test stands with subsets of the full detector readout.
//	
// Copyright Information:
//	Copyright (C) 1999	University of California, San Diego
//	Copyright (C) 2000	University of Colorado, Boulder
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include "DchCalib/DchDedxCalibList.hh"
#include "DchCalib/DchDedxCalib.hh"
#include "DchCalib/DchCalibFun.hh"
#include "boost/shared_ptr.hpp"

static const char rcsid[] = "$Id: DchDedxCalibList.cc 88 2010-01-14 12:32:57Z stroili $";

#if defined(__SUNPRO_CC) && __SUNPRO_CC <= 0x540
// FIXME: Boost 1.31 has a problem on Solaris, see if this works
#include <typeinfo>

namespace
{
   std::type_info const & ti =
     typeid( boost::checked_deleter<DchDedxCalib> );
}
#endif // Sun CC <= 5.4


DchDedxCalibList::DchDedxCalibList(
      double globalGasGain,
      std::auto_ptr<DchWireArray<double> >& gainFactors,
      std::auto_ptr<DchWireArray<DchCalibFunRCPtr> >& corrFun)
  :_globalGasGain(globalGasGain),
   _gainFactors(gainFactors),
   _corrFun(corrFun),
   _cache(*_corrFun,DchDedxCalibRCPtr((DchDedxCalib*)0)) // init _cache to all 0's, copy shape from corrFun array
{
  static bool forcePopulateCache(true);
  if (forcePopulateCache) {
          for (DchWireArray<DchDedxCalibRCPtr>::iterator i=_cache.begin();i!=_cache.end();++i)
                getDedxCalib(i.layer(),i.wire());

  }
}

DchDedxCalibList::~DchDedxCalibList()
{
}

const DchDedxCalib*
DchDedxCalibList::getDedxCalib(unsigned l,unsigned w) const
{
  if (!_cache.isValid(l,w)) return 0;
  boost::shared_ptr<DchDedxCalib>& x = _cache(l,w);
  if (x.get()==0 && (*_corrFun)(l,w)!=0 ) {
    x.reset(new DchDedxCalib(_globalGasGain,
                             (*_gainFactors)(l,w),
                             (*_corrFun)(l,w)));
  }
  return x.get();
}
