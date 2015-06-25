//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: CombinedConsistency.cc 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gautier Hamel de Monchenault
//      Alexandre Telnov, December 2007: 
//         add logLikelihood 
//
//------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "BTrk/ProbTools/CombinedConsistency.hh"

//-------------
// C Headers --
//-------------
extern "C" {
#include <assert.h>
#include <math.h>
}

//---------------
// C++ Headers --
//---------------
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BTrk/ProbTools/ConsistencySet.hh"
using std::ostream;

//----------------
// Constructors --
//----------------
CombinedConsistency::CombinedConsistency(ConsistencySet* genealogy)
  : Consistency()
  , _genealogy(genealogy)
{
  assert ( 0 != _genealogy);

  // combine consistencies
  combine();
}

CombinedConsistency::CombinedConsistency(const ConsistencySet& genealogy)
  : Consistency()
  , _genealogy(new ConsistencySet(genealogy))
{
  assert ( 0 != _genealogy);

  // combine consistencies
  combine();
}

CombinedConsistency::CombinedConsistency( const CombinedConsistency& o )
  : Consistency(o)
  , _genealogy(0)
{
  const ConsistencySet* genealogy = o.genealogy();
    
  if ( genealogy != 0 ) {
    _genealogy = new ConsistencySet( *genealogy );
  }
}

// default ctor is to be used only by subclasses that that combine the
// Consistency objects differently
CombinedConsistency::CombinedConsistency()
  : Consistency()
  , _genealogy(0)
{
}


//--------------
// Destructor --
//--------------
CombinedConsistency::~CombinedConsistency()
{
  if ( 0 != _genealogy ) {
    delete _genealogy;
  }
}

//-------------
// Methods   --
//-------------

//-------------
// Operators --
//-------------
CombinedConsistency& 
CombinedConsistency::operator=(const CombinedConsistency& rhs)
{
  if ( this != &rhs ) {
    _stat=rhs._stat;
    _value=rhs._value;
    _likelihood=rhs._likelihood;
    _logLikelihood=rhs._logLikelihood;
    _sign=rhs._sign;
    
    delete _genealogy;
    const ConsistencySet* genealogy = rhs.genealogy();
    
    if ( genealogy != 0 ) {
      _genealogy = new ConsistencySet( *genealogy );
    } else {
      _genealogy = 0;
    }
  }

  return *this;
}


//-------------
// Selectors --
//-------------
  
//-------------
// Modifiers --
//-------------

void
CombinedConsistency::combine()
{
  // copied from Consistency::combine()

  if ( 0 != _genealogy ) {

    double prodSL = 1.;
    double prodL  = 1.;
    size_t count = 0;
    double sumLogLike = 0.;

    Consistency::ConsistentStatus status = Consistency::OK;
    for(size_t i=0; i<_genealogy->nParents(); i++) {
      const Consistency* c = _genealogy->getConsistency(i);

      if (c != 0 && c->status() != Consistency::noMeasure ) {
	prodSL *= c->significanceLevel();
	prodL  *= c->likelihood();
	sumLogLike += c->logLikelihood();

	// special cases. unPhysical status wins
	if ( c->status() == Consistency::unPhysical ) {
          status = Consistency::unPhysical;
        }

      // underFlow also does
	if ( c->status() == Consistency::underFlow && 
	     status != Consistency::unPhysical ) {
	  status = Consistency::underFlow;
	}

	count++;
      }
    }
  
    // check the status
    if( prodSL>0. && count>0) {  
      // compute the combined significance level
      // formula by LeDiberder : 
      // SL = prodSL*Sum_{i=0,N-1}( (-1)^i/i! * ( ln( prodSL ) )^i )

      double   sum(1.); // the first term of the sum is unity
      double   temp(1.);
      unsigned fact(1);
      double   lnP( log( prodSL ) );
      for( size_t i=1; i<_genealogy->nParents(); i++ ) {

	temp *= -lnP;
	fact *= i;
	sum  += temp/fact;
      }
      prodSL = prodSL*sum; // significance level
    } else if (count>0) {
      prodSL = 0.;
      // at least one significance is null
      if ( status == Consistency::OK ) {
	// all statuses were OK yet one of SL values was zero
	status = Consistency::unPhysical;
      }
    } else {
      status = Consistency::noMeasure;
    } 

    // set the data members of the base class
    _value      = prodSL; // significance level
    _likelihood = prodL;
    _logLikelihood = sumLogLike;
    _sign       = Consistency::unknown;
    setStatus(status);
  }
}

void
CombinedConsistency::print(ostream& os) const
{
  Consistency::print(os);
  if ( 0 != _genealogy ) {
    _genealogy->print(os);
  }
}
