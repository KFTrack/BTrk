//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkId.hh,v 1.10 2002/03/27 22:02:57 steinke Exp $
//
// Description:
//     Holds an ID number and a pointer to an Id manager that knows how 
// provide the next number in the sequence; designed to provided tracks 
// with unique ids.  It does _not_ own the Id manager.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//------------------------------------------------------------------------
#ifndef TRKID_HH
#define TRKID_HH

// Class interface //
class TrkId {

public:
  TrkId(long idNo);    // creates with input id number
  TrkId();              // gets id number from idMan
  TrkId(const TrkId &);                    // copies existing value
  TrkId&   operator= (const TrkId&);       // copies existing value
  bool operator<(const TrkId &) const;
  ~TrkId();

  void setNewValue(const TrkId&);    // gets next Id number and copies manager
  operator long() const {return _value;}  // automatic conversion to long
  
private:	

  int _value;
  static unsigned _maxval;
};

#endif


