//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFitStatus.cc,v 1.10 2004/08/06 06:31:41 bartoldu Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
// Revision History:
//	20020502  Michael Kelsey -- Remove size argument from _history
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkFitStatus.hh"
#include "ErrLogger/ErrLog.hh"
#include <iostream>
#include <iterator>
using std::endl;
using std::ostream;

TrkFitStatus::TrkFitStatus() :
  _fitValid(false),
  _fitCurrent(false),
  _is2d(false),
  _multScat(false)
{}

TrkFitStatus::~TrkFitStatus()
{
}

// Copy ctor
TrkFitStatus::TrkFitStatus(const TrkFitStatus& right) :
  _history(right._history),
  _fitValid(right._fitValid),
  _fitCurrent(right._fitCurrent),
  _is2d(right._is2d),
  _multScat(right._multScat)
{}

TrkFitStatus&
TrkFitStatus::operator= (const TrkFitStatus& right)
{
  if(&right != this){
    _is2d = right._is2d;
    _multScat = right._multScat;
    _fitValid = right._fitValid;
    _fitCurrent = right._fitCurrent;
    _history = right._history;
  }
  return *this;
}

void
TrkFitStatus::setValid(bool v)
{
  _fitValid = v;
  if (!v) setCurrent(false);
}

ostream& operator<<( ostream& os, const TrkFitStatus& s )
{
  return s.printStatus(os);
}

ostream&
TrkFitStatus::printStatus(ostream& os) const
{
  os << " 3-d: " << (is2d() == 0) << " ";
  return printHistory(os);
}

void
TrkFitStatus::addHistory(const TrkErrCode& status,
                         const char* modulename)
{
// append the module name to the error code message
  _history.push_back(TrkHistory(status,modulename));
}

ostream&
TrkFitStatus::printHistory(ostream& os) const
{
  //FIXME: this doesn't work because ostream isn't a std::ostream, so 
  //       std::ostream_iterator doesn't like this...
  // std::copy(_history.begin(),_history.end(),std::ostream_iterator<TrkHistory>(os,"\n"));
  for(history_iterator i=_history.begin();i!=_history.end();++i)
          os << *i << endl;
  return os;
}
