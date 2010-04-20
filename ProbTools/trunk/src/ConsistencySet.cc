//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: ConsistencySet.cc 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//      Class ConsistencySet -- see header
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Yury Kolomensky       05/02/98
//
//	Copyright (C) 1998      Caltech
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "ProbTools/ConsistencySet.hh"
#include "BbrStdUtils/String.hh"

#include <algorithm>

//-------------
// C Headers --
//-------------
extern "C" {
#include <assert.h>
}

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <iomanip>
using std::endl;
using std::ostream;
using std::setw;

//----------------
// Constructors --
//----------------
ConsistencySet::ConsistencySet()
{
}

ConsistencySet::ConsistencySet( const ConsistencySet& other )
  : _consistencyList( other._consistencyList)
  , _labelList( other._labelList )
{
}

//--------------
// Destructor --
//--------------
ConsistencySet::~ConsistencySet()
{
}

//-------------
// Methods   --
//-------------
size_t
ConsistencySet::nParents() const
{
  return _labelList.size();
}

const Consistency*        
ConsistencySet::getConsistency( size_t index ) const
{
  if( index>=_consistencyList.size() ) return 0;
  return &_consistencyList[ index ];
}

const Consistency*        
ConsistencySet::getLabelConsistency( const char* label ) const
{
  size_t index;

  if ( 0 != label ) {
    std::string d(label);
    babar::String::transformToLower(d);    
    if( !getLabelIndex_( d, index ) ) return 0;
  }

  return getConsistency( index );
}

const char* 
ConsistencySet::getLabel( size_t index ) const
{
  if( index>=_labelList.size() ) return 0;
  return _labelList[ index ].c_str();
}

double      
ConsistencySet::worstSignificance() const
{
  size_t index;
  if( !worstSignificanceIndex_( index ) ) return 0.;
  const Consistency* cons = getConsistency( index );
  assert( cons!=0 );
  return cons->significanceLevel();
}
 
const char* 
ConsistencySet::worstSignificanceLabel() const
{
  size_t index;
  if( !worstSignificanceIndex_( index ) ) return 0;
  return getLabel( index );
}

ConsistencySet* 
ConsistencySet::overlap(const ConsistencySet& other) const
{
  ConsistencySet* result = 0;

  for (size_t i=0; i<this->nParents(); i++ ) {
    if ( 0 != other.getLabelConsistency(this->getLabel(i)) ) {
      if ( 0 == result ) {
	result = new ConsistencySet;
      }

      result->add(this->getLabel(i), *(this->getConsistency(i)) );
    }
  }

  return result;
}



//-------------
// Operators --
//-------------
ConsistencySet& 
ConsistencySet::operator=(const ConsistencySet& rhs)
{
  if ( this != &rhs ) {
    _consistencyList=rhs._consistencyList;
    _labelList = rhs._labelList;
  }

  return *this;
}

bool 
ConsistencySet::operator==(const ConsistencySet& rhs) const {
  bool retval = false;

  // just in case
  if ( this != &rhs ) {
    
    if ( this->nParents() == rhs.nParents() ) {
      
      retval = true;
      for (size_t i=0;i<nParents();i++) {
	// order is not guaranteed to be the same even for equal objects
	if ( 0 == rhs.getLabelConsistency(this->getLabel(i)) ) {
	  retval = false;
	  return retval;
	}
      }
    }
  }

  return retval;
}

bool 
ConsistencySet::operator!=(const ConsistencySet& rhs) const {
  return !(*this == rhs);
}

//-------------
// Selectors --
//-------------
  
//-------------
// Modifiers --
//-------------

void 
ConsistencySet::reset()
{
  _consistencyList.clear();
  _labelList.clear();
}

bool
ConsistencySet::add( const char* label, const Consistency& c )
{
  bool retval = false;

  // some idiot-proofing

  if ( 0 != label ) {

    std::string d(label);
    babar::String::transformToLower(d);

    // test for overlaps
    size_t index;

    if( ! getLabelIndex_( d, index ) ) {
      _labelList.push_back( d );
      _consistencyList.push_back( c );
      retval = true;
    }
  }

  return retval;
}

bool 
ConsistencySet::combine(const ConsistencySet& other) {
  bool retval = true;
  for (size_t i=0;i<other.nParents();i++) {
    retval = retval && add(other.getLabel(i),*(other.getConsistency(i)));
  }

  return retval;
}


//
//  Internal protected functions
//

bool
ConsistencySet::getLabelIndex_( const std::string& dLower, 
				size_t& index ) const
{
  bool retval = false;

  std::vector<std::string>::const_iterator iter = _labelList.begin();
  index = 0;
  while (iter != _labelList.end()) {
    if (*iter == dLower) {
      retval = true;
      break;
    }
      index++;
      ++iter;
  }
    
  return retval;
}

bool   
ConsistencySet::worstSignificanceIndex_( size_t& index ) const
{
  double worst( 999. );
  const Consistency* cons;
  size_t i=0;
  while( cons=getConsistency(i) )
    {
      double value=cons->significanceLevel();
      if( value<worst )
	{
	  index=i;
	  worst=value;
	}
      i++;
    }
  return worst<=1.;
}


void
ConsistencySet::print(ostream& os) const
{
  os << nParents() << " parents:" << endl;
  for (size_t i=0;i<nParents();i++) {
    os << setw(3) << i;
    os << setw(10) << _labelList[i] << " : ";
    _consistencyList[i].print(os);
  }
}
  
