#ifndef DCHDBIOPARSER_HH
#define DCHDBIOPARSER_HH

//--------------------------------------------------------------------------
//
// Environment:
//      This software was developed for the BaBar collaboration.  If you
//      use all or part of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998      <Institution>
//
//------------------------------------------------------------------------

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//----------------------
// Base Class Headers --
//----------------------
#include "DchGeomBase/DchGFileParser.hh"

//		---------------------
// 		-- Class Interface --
//		---------------------


/**
 *  C++ source file code DchDbioParser. The first sentence is a brief summary of 
 *  what the class is for. It is followed by more detailed information
 *  about how to use the class. This doc comment must immediately preced the 
 *  class definition.
 *
 *  Additional paragraphs with more details may follow; separate paragraphs
 *  with a blank line. The last paragraph before the tags (preceded by @) 
 *  should be the identification and copyright, as below.
 *
 *  Please note that KDOC comments must start with a forward slash
 *  followed by TWO asterisks. Interface members should be documented
 *  with KDOC comments, as should be protected members that may be of interest
 *  to those deriving from your class. Private implementation should
 *  be commented with C++-style // (double forward slash) comments.
 *
 *  This software was developed for the BaBar collaboration.  If you
 *  use all or part of it, please give an appropriate acknowledgement.
 *
 *  Copyright (C) 1998 [your institution]
 *
 *  @see DchDbioParserDchDbioParser
 *
 *  @version $Id: DchDbioParser.hh 92 2010-01-14 12:38:30Z stroili $ 
 *
 *  @author (Author1)		(originator/contributor etc.);
 *  @author (Author2)		(originator/contributor etc.)
 */

class DchDbioParser : public DchGFileParser {

//--------------------
// Declarations     --
//--------------------

    // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DchDbioParser( const char * );

  // Destructor
  virtual ~DchDbioParser( );
  
  // Operators
    
  //    virtual int operator==( const DchDbioParser& ) const;
  //            int operator!=( const DchDbioParser& ) const;
  
  /**
   *  This is a member function with an non-obvious purpose or interface.
   *  Use a KDOC comment to describe its use. Use the @see tag to add
   *  references to other members (of this or another class) or classes.
   *
   *  Use descriptive rather than imperative language (i.e.
   *  "Performs an important operation" rather than "Perform an important
   *  operation").
   *
   *  @param myParameter  a parameter supplied to the function
   *  @return             the return code, zero for success 
   *
   *  @see DchDbioParserDchDbioParser#myFunction
   */
//   int myFunction( const DchDbioParserDchDbioParser& myParameter );
  
  // Selectors (const)
  //  get Dch volumes
  bool getDchVol( const std::string& volName, 
		  std::vector<double>& par,  
		  std::vector<double>& position );
  //  get Dch sense wires
  bool getDchSWire( int lay, std::vector<double>& par);  
  //  get Dch field wires
  bool getDchFWire( const std::string& type, int lay, 
		    std::vector<double>& par);  

  // Modifiers
  
//   enum sense {cells=0, Srad, Sphi, Stwist, Sstereo, Sdphi, Ssag, Sdiam, Svolt};
//   enum iofield {nwires=0, Fradius, Fphi, Ftwist, Fstereo, Fdphi, Fsag, Fdiam, 
// 		Fvolt};
//   enum lfield {LFrad1=0, LFrad2, LFphi, LFtwist, LFstereo1, LFstereo2, LFdphi, 
// 	       LFsag, LFdiam, LFvolt};
  
private:
  
  // Data members
//   ifstream _inFile;             // input file stream
  std::streampos _begin;
  std::string _filename;
  // Friends
  
  
  // Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.
  DchDbioParser( const DchDbioParser& );                // Copy Constructor
  DchDbioParser& operator= ( const DchDbioParser& );    // Assignment op
  
  bool seek( const std::string& );
  bool seekTemplate( const std::string& );
  bool locate( const std::string& templ, const std::string& name);
  bool check( const std::string&, const std::string& name );

};

#endif // DCHDBIOPARSER_HH
