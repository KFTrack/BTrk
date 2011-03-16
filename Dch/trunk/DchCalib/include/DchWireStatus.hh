//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchWireStatus.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchWireStatus.
//
//      Simplistic class which provides an interface (and nothing 
//      more than an interface) to the information 
//      derived by the rolling calibration classes in OPR.
//      For documentation on how these numbers are derived,
//      check DchOpr/DchOprEffTrk and DchOpr/DchChanNoise and
//      http://www.slac.stanford.edu/~jalbert/talk4_16_99.ps
//      http://www.slac.stanford.edu/~jalbert/talkjune4.pdf
//      http://www.slac.stanford.edu/~jalbert/effplots6_6_99.ps
//
//      The 'isNoisy' and 'isDead' bools map to the infomation 
//      in the Dch 'CalStatusChan'.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gerhard Raven           11/01/99
//
// Copyright Information:
//      Copyright (C) 1999      University of California, San Diego
//
//------------------------------------------------------------------------

#ifndef DCHWIRESTATUS_HH
#define DCHWIRESTATUS_HH


class DchWireStatus 
{
public:
        ~DchWireStatus();

        unsigned layer() const { return _layer; }
        unsigned wire() const { return _wire; }

        bool isOK() const { return !isNoisy() && !isDead();}
        bool isNoisy() const { return _noisy;}
        bool isDead() const  { return _dead;}
        bool hasBadCharge() const { return _badCharge; }

        double efficiency() const;
        double efficiencyErr() const;
        double noise() const;
        double noiseErr() const;


private:
        friend class DchWireStatusList;
        friend class DchWireStatusListProxy;
        friend class DchCdbWireStatusListProxy;

        DchWireStatus(unsigned layer, unsigned wire,
                      bool dead, bool noisy, bool hasBadCharge,
                      double eff=-1,double effErr=-1,
                      double noise=-1,double noiseErr=-1);

        unsigned _layer,_wire;
        double _eff,_effErr;
        double _noise,_noiseErr;
        bool _noisy,_dead, _badCharge;

};

#endif // DCHWIRESTATUS_HH
