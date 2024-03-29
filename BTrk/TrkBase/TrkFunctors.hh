//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFunctors.hh,v 1.3 2004/09/30 16:59:41 hulsberg Exp $
//
// Description:
//      general, simple tracking predicates
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gerhard Raven
//------------------------------------------------------------------------
#ifndef TRKFUNCTORS_HH
#define TRKFUNCTORS_HH
#include <functional>

#include "BTrk/TrkBase/TrkHit.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
class TrkHitUpdater;
class TrkRep;

namespace TrkBase { namespace Functors {


        template <class T> struct takeAddress {
                T* operator()(T& t) { return &t; }
        };


        // Functors which are only accessible by inheriting from TrkHitUpdater
        // FIXME: maybe these should live in the TrkHitUpdater namespace instead???

        class updateMeasurement {
        public:
                TrkErrCode operator()(TrkHit& h) const
                { return h.updateMeasurement(_t); }
        private:
                // Only TrkHitUpdater can create one of these...
                friend class ::TrkHitUpdater;
                updateMeasurement(const TrkDifTraj* traj=0) : _t(traj) {}
                const TrkDifTraj *_t;
        };

        class setActive {
        public:
                TrkHit* operator()(TrkHit& h) const
                { return h.setActive(_a); }
                TrkHit* operator()(TrkHit*& h) const
                { return h->setActive(_a); }
        private:
                // Only TrkHitUpdater can create one of these...
                friend class ::TrkHitUpdater;
                setActive(bool active) : _a(active){}
                bool _a;
        };

        class setParent {
        public:
                TrkHit* operator()(TrkHit& h) const
                { return h.setParent(_p); }
        private:
                // Only TrkHitUpdater can create one of these...
                friend class ::TrkHitUpdater;
                setParent(TrkRep* parent) : _p(parent){}
                TrkRep* _p;
        };
}}
#endif
