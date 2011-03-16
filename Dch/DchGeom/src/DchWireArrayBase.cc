#include "BaBar/BaBar.hh"
#include "DchGeom/DchWireArrayBase.hh"
#include "DchGeom/DchDetector.hh"
#include "DchGeom/DchLayout.hh"

#include <functional>
#include <algorithm>


DchWireArrayBase::DchWireArrayBase(const DchDetector& det)
                : _first(det.nLayer(),0),    //FIXME: at some point DchDetector should tell us 
                  _offset(det.nLayer()+1,0) // what the first wire in each layer actually is...
{
        for (unsigned l=1;l<_offset.size();++l)
                _offset[l] = _offset[l-1]+det.nWires(l);
}

DchWireArrayBase::DchWireArrayBase(const DchLayout& det)
                : _first(det.nLayers(),0),    //FIXME: at some point DchLayout should tell us 
                  _offset(det.nLayers()+1,0) // what the first wire in each layer actually is...
{
        for (unsigned l=1;l<_offset.size();++l)
                _offset[l] = _offset[l-1]+det.nWires(l);
}

void
DchWireArrayIteratorBase::inc(const DchWireArrayBase& p, int i) {
        assert(i>0);
        if(atEnd(p)) return;
        _w+=i;
        while(_w>p.lastWire(_l)&&_l<p.lastLayer()) {
                _w-=p.lastWire(_l);
                _w+=p.firstWire(++_l)-1;
        }
        if(_l==p.lastLayer()) _w=std::min(_w,p.firstWire(_l));
}

void
DchWireArrayIteratorBase::dec(const DchWireArrayBase& p, int i) {
        assert(i>0);
        if (atBegin(p)) return;
        _w-=i;
        while(_w<p.firstWire(_l)&&_l>p.firstLayer()) {
                _w+=p.firstWire(_l);
                _w-=p.lastWire(--_l)+1;
        }
        if(_l==p.firstLayer()) _w=std::max(_w,p.firstWire(_l));
}
