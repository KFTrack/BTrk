//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkExtInterface.hh,v 1.9 2000/12/05 19:59:31 brownd Exp $
//
// Description:
//     Abstract base class; derived classes enable TrkRecoTrk to give 
// users access to representation-specific functions without having 
// the functions in the RecoTrk interface.  Users create an object of 
// of one of the derived TrkExtInterface classes and pass it to 
// TrkRecoTrk::attach(), which will either accept it or reject it, 
// depending on whether the TrkRep inside recognizes it.  Recognition 
// is implemented by a simple equality test on the Rep's key.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#ifndef TRKEXTINTERFACE_HH
#define TRKEXTINTERFACE_HH

class TrkRep;
class IfdKey;

// Class interface //
class TrkExtInterface {

public:
  TrkExtInterface();
  virtual ~TrkExtInterface();
  virtual bool attach(TrkRep*);
  virtual bool attach(const TrkRep*);
  bool isAttached() const                                 {return _myRep != 0;}
  bool nonConstAttachment() const { return _nonconst; }
protected:	
  TrkRep* myRep();
  const TrkRep* myConstRep() const;
  virtual const IfdKey& myKey() const = 0;
  void setRep(const TrkRep*);
  void setRep(TrkRep*);
private:	
  TrkRep* _myRep;
  bool _nonconst; // keep track of whether we attached const or not
// Preempt 
  TrkExtInterface&   operator= (const TrkExtInterface&);
  TrkExtInterface(const TrkExtInterface &);
};

#endif







