//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFitStatus.hh,v 1.11 2004/08/06 06:31:41 bartoldu Exp $
//
// Description:
//     Describes the status of a track fit.  Designed to be inherited by 
// TrkReps, and can be used to present one facet of them to Trk users.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#ifndef TRKFITSTATUS_HH
#define TRKFITSTATUS_HH

#include "BTrk/TrkBase/TrkHistory.hh"
#include <vector>
#include <utility>
#include <functional>
#include <iostream>

// Class interface //
class TrkFitStatus {
  typedef std::vector<TrkHistory>::const_iterator history_iterator;
  typedef std::vector<TrkHistory>::const_reverse_iterator history_riterator;

public:
  bool is2d()       const      {return _is2d;}
  bool fitCurrent() const      {return _fitCurrent;}
  bool fitValid()   const      {return _fitValid;}
  bool multScat()   const      {return _multScat;}      
  std::ostream& printStatus(std::ostream& os=std::cout) const;

  void setValid(bool v);
  void setCurrent(bool c) { _fitCurrent = c; }
  void set2d(bool d) { _is2d = d; };
  void setMultScat(bool m) { _multScat = m; };
// here's the error code returned by the most recent fit
  const TrkErrCode& fitStatus() const { return _history.back().status(); }
// a complete history of how the track was created/modified
//
  history_iterator beginHistory() const { return _history.begin(); }
  history_iterator endHistory() const   { return _history.end(); }
  history_riterator reverseBeginHistory() const { return _history.rbegin(); }
  history_riterator reverseEndHistory() const   { return _history.rend(); }
  std::pair<history_iterator,history_iterator> history() const {
          return std::pair<history_iterator,history_iterator>(beginHistory(),endHistory());
  }
  const std::vector<TrkHistory>& historyVector() const { return _history; }
// add to the tracks history.  Provide the module name, please!
  virtual void addHistory(const TrkErrCode& status,const char* modulename);
  template <class T>
  void addHistory(T begin, T end) { _history.insert(_history.end(),begin,end); }
  template <class T>
  void addHistory(std::pair<T,T> p) { addHistory(p.first,p.second); }
  std::ostream& printHistory(std::ostream& os=std::cout) const;
protected:
  virtual ~TrkFitStatus();
  TrkFitStatus();
  TrkFitStatus&   operator= (const TrkFitStatus&);
  TrkFitStatus(const TrkFitStatus &);
private:
  std::vector<TrkHistory> _history;
  bool _fitValid;
  bool _fitCurrent;
  bool _is2d;
  bool _multScat;
};

std::ostream& operator<<( std::ostream& os, const TrkFitStatus& s );


#endif

