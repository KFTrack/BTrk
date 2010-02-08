// $Id: DetElemPrinter.cc,v 1.11 2007/01/27 07:49:03 wiredces Exp $

#include "BaBar/BaBar.hh"

#include "DetectorModel/DetElemPrinter.hh"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>

#include "DetectorModel/DetType.hh"
#include "DetectorModel/DetElem.hh"
#include "DetectorModel/DetElemPointIterator.hh"
#include "DetectorModel/DetSet.hh"
using std::endl;
using std::ios;
using std::ostream;
using std::setiosflags;
using std::setprecision;
using std::setw;

void
DetElemPrinter::print( ostream& o ) const
{
  
  // get the list of all Dirc elements
  DetElemList theElemList;
  // drc.drcSet().listAllElements( theElemList );
  _set.listAllElements( theElemList );
  
  // create the hash dictionary
  std::map< const DetType*, DetElemList* > theHash;
  //&detElemPrinterTypeHash );// This used to be an argument to theHash.(?)

  // Loop over Dirc elements
  for ( DetElemList::const_iterator theElemIterator = theElemList.begin();
	theElemIterator != theElemList.end();
	++theElemIterator )
    {
      DetElem* elem = *theElemIterator;
      const DetType* type = elem->detectorType();
      if (theHash[type] == NULL) {
	theHash[type] = new DetElemList();
      }
      (theHash[type])->push_back( elem ) ;
    }

  for ( std::map< const DetType*, DetElemList* >::const_iterator theHashIter =
	  theHash.begin();
	theHashIter != theHash.end();
	++theHashIter )
    {
      DetElemList* elemList = (*theHashIter).second;
      const DetType* type = (*theHashIter).first;
      if( type!=0 ) 
	{
	  o << endl << "******* NEW TYPE *******" << endl;
	  std::string s( type->typeName() );
	  removeBlanks( s );	
	  o << s << "  " << type->typeNumber() << "  ";
	}
      else
	{
	  continue;
	  //	  o << endl << "******* NULL TYPE ******" << endl;
	}
      o << "nElem " << elemList->size() << endl;
      if( type!=0 )
	{
	  for ( DetElemList::const_iterator elemIter = elemList->begin();
		elemIter != elemList->end();
		++elemIter )
 	    {
	      const DetElem* elem = *elemIter;
	      std::string s( elem->elementName() );
	      removeBlanks( s );
	      o << s << "  " << elem->elementNumber() << endl;

	      DetElemPointIterator iter( *elem );
	      HepPoint aVertex;
	      std::vector< HepPoint > vertexList;
	      bool move=false;

	      Action firstPoint = iter.next( aVertex );
	      vertexList.push_back( aVertex );

	      Action whatNext = iter.next( aVertex );
	      if ( whatNext == Continue ) 
		{
		  while(1) 
		    {
		      if( whatNext==CloseLine || whatNext==CloseShape ) 
			move= true;
		      if( whatNext == CloseShape ) break;
		      if( whatNext == CloseLine  ) move=true;
		      if (move)  
			{
			  printFace( vertexList, o );
			  vertexList.clear();
			}
		      vertexList.push_back(aVertex);
		      move=false;
		      whatNext = iter.next( aVertex );
		    }
		}
	      // if vertexList not empty, draw corresponding FaceSet
	      if( vertexList.size()>0 ) 
		{
		  printFace( vertexList, o );
		  vertexList.clear();
		}
	    }
	}
    }
}

void
DetElemPrinter::printFace( const std::vector< HepPoint > & vertexList, ostream& o ) const
{
  o << setiosflags(ios::fixed) << setw(4);
  o <<  vertexList.size() << "  ";
  for( size_t i=0; i<vertexList.size(); i++ )
    {
      HepPoint point = vertexList[i];
      o << setiosflags(ios::fixed) << setw(8) << setprecision(2);
      o << point.x() << " ";
      o << point.y() << " ";
      o << point.z() << " ";
      o << "     ";
    } 
  o << endl;
}

void
DetElemPrinter::removeBlanks( std::string & aString ) const
{
  aString.erase( aString.find_last_not_of(" ") + 1 ); // trailing
  aString.erase( 0, aString.find_first_not_of(" ") ); // leading
  size_t i;

  while( 1 )
    {
      i=aString.find_first_of(" ");

      if( i == std::string::npos )
        {
          break;
        }
      aString.replace( i, 1, "_" );
    }
}

unsigned 
detElemPrinterTypeHash( const DetType* const & p )
{
  unsigned long i = (unsigned long)p;
  return  ( i >> 5) & 0x03ff;
}
