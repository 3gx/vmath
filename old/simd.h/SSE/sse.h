#pragma once

#ifndef __SSE__
#error "Please enable SSE support in g++"
#endif

#define __out
#include <cassert>

namespace SIMD
{
  typedef float  _v4sf  __attribute__((vector_size(16)));
  typedef double _v2df  __attribute__((vector_size(16)));

#include "sse_fp32.h"
#include "sse_fp64.h"
};

inline SIMD::dvec  LDA(const double &x) {return SIMD::load<true >(x);}
inline SIMD::dvec  LDU(const double &x) {return SIMD::load<false>(x);}
inline SIMD::dvec& STA(double &x) {return *(SIMD::dvec*)&x;}
inline SIMD::dvec STU(double &x, const SIMD::dvec y) {SIMD::store<false>(x, y); return y;}

inline SIMD::fvec  LDA(const float &x) {return SIMD::load<true >(x);}
inline SIMD::fvec  LDU(const float &x) {return SIMD::load<false>(x);}
inline SIMD::fvec& STA(float &x) {return *(SIMD::fvec*)&x;}
inline SIMD::fvec STU(float &x, const SIMD::fvec y) {SIMD::store<false>(x, y); return y;}

