//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkIdManager.hh,v 1.5 2003/04/04 00:33:40 brownd Exp $
//
// Description:
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------
#ifndef TRKIDMANAGER_HH
#define TRKIDMANAGER_HH
#include "AbsEvent/AbsEvtObj.hh"

// Class interface //
class TrkIdManager : public AbsEvtObj {

public:
  TrkIdManager();
  virtual ~TrkIdManager();

  virtual long nextId() = 0;       // Get next id number and update manager
  virtual long lastId() const = 0;
  virtual void setMax(long maxid) = 0;

private:	
  // Preempt 
  TrkIdManager&   operator= (const TrkIdManager&);
  TrkIdManager(const TrkIdManager &);
};

#endif
