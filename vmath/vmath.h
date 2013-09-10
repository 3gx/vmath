#pragma once
#include <cassert>

/* function defintions */

/************************/
/****** INTEGER *********/
/************************/

template<typename vint> static inline vint vspan(const int);

/***** LOAD/STORE *******/

/* vector load */
template<typename vint> static inline vint vload(const int*);
/* vector stores */
template<typename vint> static inline vint vstore(int*, const vint);
/* gather load */
template<typename vint, typename vint_addr> static inline vint vgather(const int*, const vint_addr);
/* scatter store */
template<typename vint, typename vint_addr> static inline vint vscatter(int*, const vint_addr, const vint_addr);

/******** MATH ***********/

/* min */
/* max */
/* abs */

/************************/
/******* DOUBLE *********/
/************************/

/***** LOAD/STORE *******/

/* vector load */
template<typename vreal> static inline vreal vload(const double*);
/* vector stores */
template<typename vreal> static inline vreal vstore(double*, const vreal);
/* gather load */
template<typename vreal, typename vint> static inline vreal vgather(const double*, const vint);
/* scatter store */
template<typename vreal, typename vint> static inline vreal vscatter(double*, const vint, const vreal);

/******** MATH ***********/

/* min */
/* max */
/* abs */
/* pow */
/* exp */
/* log */
/* log10 */
/* sqrt */
/* invsqrt  */
/* cube root  */
/* inverse cube root  */
/* cos */
/* sin */
/* square */
template<typename vreal> 
static inline vreal __square(const vreal x)  { return x*x; }
/* integer pow */
template<int POW,typename vreal> static inline vreal __ipow(const vreal x)
{
  static_assert(POW >= 1 && POW <= 6, "1 <= POW <= 6 failed! Please check your code");
  switch(POW)
  {
    case 1: return x;
    case 2: return __square(x);
    case 3: return __square(x)*x;
    case 4: return __square(__ipow<2>(x));
    case 5: return __square(__ipow<2>(x))*x;
    case 6: return __square(__ipow<3>(x));
    default: assert(0);
  };
}

/********* CONDITIONALS ***********/

template<typename vmask> static inline vmask vtrue ();
template<typename vmask> static inline vmask vfalse();
template<typename vreal> struct Mask;
template<typename T> static inline bool all (const T condition) { return Mask<T>(condition).all (); }
template<typename T> static inline bool any (const T condition) { return Mask<T>(condition).any (); }
template<typename T> static inline bool none(const T condition) { return Mask<T>(condition).none(); }
  
  template<typename T, typename _If, typename _ElseIf, typename _Else> 
static inline void cifelseif(const Mask<T> &mask, const T &condition_if, const _If &If, const T &condition_elseif, const _ElseIf &ElseIf, const _Else &Else) __attribute__((always_inline))
{
  const Mask<T> maskif    ( condition_if                       && (T)mask);
  const Mask<T> maskelseif(!condition_if &&  condition_elseif  && (T)mask);
  const Mask<T> maskelse  (!condition_if && !condition_elseif  && (T)mask);

  if     (maskif    .all ()) { If    (maskif    );                     }
  else if (maskelseif.all ()) { ElseIf(maskelseif);                     }
  else if (maskelse  .all ()) { Else  (maskelse  );                     }
  else if (maskif    .none()) { ElseIf(maskelseif); Else  (maskelse  ); }
  else if (maskelseif.none()) { If    (maskif    ); Else  (maskelse  ); }
  else if (maskelse  .none()) { If    (maskif    ); ElseIf(maskelseif); }
  else 
  {
    If    (maskif    );
    Else  (maskelse  );
    ElseIf(maskelseif);
  }
}
  template<typename T, typename _If, typename _Else> 
static inline void cif(const Mask<T> &mask, const T &condition, const _If &If, const _Else &Else) __attribute__((always_inline))
{
  const Mask<T> maskif  ( condition && (T)mask);
  const Mask<T> maskelse(!condition && (T)mask);
  if      (maskif  .all()) { If(maskif);                 }
  else if (maskelse.all()) {             Else(maskelse); }
  else                     { If(maskif); Else(maskelse); }
}
  template<typename T, typename _If>
static inline void cif(const Mask<T> &mask, const T &condition, const _If &If) __attribute__((always_inline))
{
  const Mask<T> maskif(condition && (T)mask);
  if (maskif.any()) If(maskif);
}
 
  template<typename T, typename _If, typename _ElseIf, typename _Else> 
static inline void cifelseif(const T &condition_if, const _If &If, const T &condition_elseif, const _ElseIf &ElseIf, const _Else &Else) __attribute__((always_inline))
{
  cifelseif(Mask<T>(vtrue<T>()), condition_if, If, condition_elseif, ElseIf, Else);
}
  template<typename T, typename _If, typename _Else> 
static inline void cif(const T &condition, const _If &If, const _Else &Else) __attribute__((always_inline))
{
  cif(Mask<T>(vtrue<T>()), condition, If, Else);
}
  template<typename T, typename _If>
static inline void cif(const T &condition, const _If &If) __attribute__((always_inline))
{
  cif(Mask<T>(vtrue<T>()), condition, If);
}

#define vassert(condition,mask) assert(all((condition) || !(mask)))
#define IF(_mask_) [&](const Mask<vmask> &_mask_) __attribute__((always_inline))
#define ELSEIF(_mask_) [&](const Mask<vmask> &_mask_) __attribute__((always_inline))
#define ELSE(_mask_) [&](const Mask<vmask> &_mask_) __attribute__((always_inline))

template<int Step>
struct Vectorize
{
  template<typename Lambda>
    static void loop(const int nbeg, const int nend, Lambda func)
    {
      const int n = nend - nbeg;
      const int n1 = (n/Step)*Step;
      for (int ib = 0; ib < n1; ib += Step)
      {
        _Pragma("simd")
          for (int lane = 0; lane < Step; lane++)
            func(nbeg + ib+lane,lane);
      }

      for (int i = n1; i < n; i++)
        func(i,0);
    }
};


#ifdef __AVX512__
#error "AVX512 is not supported"
#include "avx512.h"
#elif defined __MIC__
#error "MIC is not supported"
#include "mic.h"
#elif defined __AVX__
#warning "AVX is enabled"
#include "avx.h"
#elif defined __SSE__
#warning "SSE is enabled"
#include "sse.h"
#else
#warning "No Vector ISA enabled, code may not compile..."
#define NOVECTORISA
#endif
#include "scalar.h"

#define VECTORLOOPFUNCTION(vreal,vint,vmask,index) \
    template<typename vreal, typename vint, typename vmask>\
      inline void operator() (const int index,\
          const vreal __fdum__, const vint __idum__, const vmask __bdum__) __attribute__((always_inline))

struct foreach
{
  template<typename Lambda>
    static void loop(const int beg, const int end, Lambda &func)
    {
#ifndef NOVECTORISA
      const int n    = end - beg;
      const int nvec = vdouble::round(n);
      for (int i = 0; i < nvec; i += vdouble::SIZE)
        func(i, vdouble(0.0), vinteger(0), vboolean(0.0));
#else
      const int nvec = 0;
#endif
      for (int i = nvec; i < n; i++)
        func(i, 0.0, 0, false);
    }
};
