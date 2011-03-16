#ifndef DCHGFILEPARSER_HH
#define DCHGFILEPARSER_HH

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

//-------------
// C Headers --
//-------------
extern "C" {
}

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

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------


/**
 *  C++ source file code DchGFileParser. The first sentence is a brief summary of 
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
 *  @see DchGFileParserDchGFileParser
 *
 *  @version $Id: DchGFileParser.hh 92 2010-01-14 12:38:30Z stroili $ 
 *
 *  @author (Author1)		(originator/contributor etc.);
 *  @author (Author2)		(originator/contributor etc.)
 */

class DchGFileParser {

//--------------------
// Declarations     --
//--------------------

    // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DchGFileParser( const char * );

  // Destructor
  virtual ~DchGFileParser( );
  
  // Operators
    
  //    virtual int operator==( const DchGFileParser& ) const;
  //            int operator!=( const DchGFileParser& ) const;
  
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
   *  @see DchGFileParserDchGFileParser#myFunction
   */
//   int myFunction( const DchGFileParserDchGFileParser& myParameter );
  
  // Selectors (const)
  //  get Dch volumes
  virtual bool getDchVol( const std::string& volName, 
			  std::vector<double>& par,  
			  std::vector<double>& position ) = 0;
  //  get Dch sense wires
  virtual bool getDchSWire( int lay, std::vector<double>& par) = 0;  
  //  get Dch field wires
  virtual bool getDchFWire( const std::string& type, int lay, 
			    std::vector<double>& par) = 0;  

  // Modifiers
  
  enum sense {cells=0, Srad, Sphi, Stwist, Sstereo, Sdphi, Ssag, Sdiam, Svolt};
  enum iofield {nwires=0, Fradius, Fphi, Ftwist, Fstereo, Fdphi, Fsag, Fdiam, 
		Fvolt};
  enum lfield {LFrad1=0, LFrad2, LFphi, LFtwist, LFstereo1, LFstereo2, LFdphi, 
	       LFsag, LFdiam, LFvolt};

protected:
  
  // Helper functions
  std::ifstream _inFile;             // input file stream
  void fileRewind( void ) { _inFile.clear(); _inFile.seekg(0); }

private:
  
  // Data members
  std::streampos _begin;
  std::string _filename;
  // Friends
  
  
  // Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.
  DchGFileParser( const DchGFileParser& );                // Copy Constructor
  DchGFileParser& operator= ( const DchGFileParser& );    // Assignment op
  
};

#endif // DCHGFILEPARSER_HH
