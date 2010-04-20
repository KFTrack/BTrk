//--------------------------------------------------------------------------
// File and Version Information:
//
// Description:
//
// Environment:
//
// Author List:
//       Gerhard Raven
//
// Copyright Information:
//
//------------------------------------------------------------------------

#ifndef BBRPAIRUTILS_HH
#define BBRPAIRUTILS_HH

#include <functional>

namespace babar {
  namespace Pair {

        // select1st and select2nd are extensions: they are not part of the standard.
        template <class _Pair>
        struct _Select1st : public std::unary_function<_Pair, typename _Pair::first_type> {
          const typename _Pair::first_type& operator()(const _Pair& __x) const {
            return __x.first;
          }
        };

        template <class _Pair>
        struct _Select2nd : public std::unary_function<_Pair, typename _Pair::second_type>
        {
          const typename _Pair::second_type& operator()(const _Pair& __x) const {
            return __x.second;
          }
        };

        template <class _Pair> struct select1st : public _Select1st<_Pair> {};
        template <class _Pair> struct select2nd : public _Select2nd<_Pair> {};

        // project1st and project2nd are extensions: they are not part of the standard
        template <class _Arg1, class _Arg2>
        struct _Project1st : public std::binary_function<_Arg1, _Arg2, _Arg1> {
          _Arg1 operator()(const _Arg1& __x, const _Arg2&) const { return __x; }
        };

        template <class _Arg1, class _Arg2>
        struct _Project2nd : public std::binary_function<_Arg1, _Arg2, _Arg2> {
          _Arg2 operator()(const _Arg1&, const _Arg2& __y) const { return __y; }
        };

        template <class _Arg1, class _Arg2> struct project1st : public _Project1st<_Arg1, _Arg2> {};
        template <class _Arg1, class _Arg2> struct project2nd : public _Project2nd<_Arg1, _Arg2> {};

  }
}
#endif
