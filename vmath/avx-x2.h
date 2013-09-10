#include <immintrin.h>
#include <iostream>

#include "vmath.h"

/* avx implementation of vmath */

#ifndef __AVX__
#error "AVX support must be enabled!"
#endif




struct int8
{
  enum {SIZE = 8};
  __m128i i1, i2;
  
  inline int8() {}
  inline int8(const __m128i v1, const __m128i v2) : i1(v1), i2(v2) {}
  inline int8(const int  s) : i1(_mm_set1_epi32(s)), i2(_mm_set1_epi32(s)) {}

  const inline int8& operator+=(const int8 a) {i1 = _mm_add_epi32(i1,a.i1); i2 = _mm_add_epi32(i2,a.i2); return *this;}
  const inline int8& operator-=(const int8 a) {i1 = _mm_sub_epi32(i1,a.i1); i2 = _mm_sub_epi32(i2,a.i2); return *this;}
  const inline int8& operator*=(const int8 a) {i1 = _mm_mullo_epi32(i1,a.i1); i2 = _mm_mullo_epi32(i2,a.i2); return *this;}
  const inline int8& operator/=(const int8 a) {i1 = _mm_div_epi32(i1,a.i1); i2 = _mm_div_epi32(i2,a.i2); return *this;}

  friend inline int8 operator + (int8 a, const int8 b) { return a+=b; };
  friend inline int8 operator - (int8 a, const int8 b) { return a-=b; };
  friend inline int8 operator * (int8 a, const int8 b) { return a*=b; };
  friend inline int8 operator / (int8 a, const int8 b) { return a/=b; };

};

/***** LOAD/STORE *******/

/******** MATH ***********/

static inline int8 __max(const int8 a, const int8 b) { return int8(_mm_max_epi32(a.i1,b.i1), _mm_max_epi32(a.i2,b.i2)); }
static inline int8 __min(const int8 a, const int8 b) { return int8(_mm_min_epi32(a.i1,b.i1), _mm_min_epi32(a.i2,b.i2)); }
static inline int8 __abs(const int8 a              ) { return int8(_mm_abs_epi32(a.i1),_mm_abs_epi32(a.i2));}

/************************/
/******* DOUBLE *********/
/************************/

static inline __m256d get_mask(const unsigned long long mask) { return(__m256d)_mm256_set1_epi64x(mask); }
static inline __m256d abs_mask() { return get_mask(~(1ULL<<63)); }
static inline __m256d full_mask() { return get_mask(-1); }
static inline __m256d empty_mask() { return get_mask(0); }

/* vector double */
struct double8
{
  enum {SIZE = 8};
  inline static int round(const int n) { return n & (-SIZE); }


  __m256d x1, x2;
  inline double8() {}
  inline double8(const __m256d v1, const __m256d v2) : x1(v1), x2(v2) {}
  inline double8(const double  s) : x1(_mm256_set1_pd(s)), x2(_mm256_set1_pd(s)) {}
  inline double8(const double *data, const int idx[8]) 
  {
    x1 = _mm256_set_pd(data[idx[3]], data[idx[2]], data[idx[1]], data[idx[0]]);
    x2 = _mm256_set_pd(data[idx[7]], data[idx[6]], data[idx[5]], data[idx[4]]);
  }
  inline double8(const double *data, const int8 _idx) 
  {
    union { __m128i v[2]; int s[8]; } idx;
    idx.v[0] = _idx.i1;
    idx.v[1] = _idx.i2;
    *this = double8(data, idx.s);
  }

  inline double8(const double *data) : x1(_mm256_loadu_pd(data)), x2(_mm256_loadu_pd(data+4)) {}

  inline operator int8() const {return int8(_mm256_cvttpd_epi32(x1), _mm256_cvttpd_epi32(x2));}

  /* basic arithmetics */
  const inline double8& operator =(const double8 a) {x1  = a.x1; x2 = a.x2; return *this;}
  const inline double8& operator =(const double  a) {x1 = x2 =  _mm256_set1_pd(a); return *this;}
  const inline double8& operator+=(const double8 a) {x1 = _mm256_add_pd(x1,a.x1); x2 = _mm256_add_pd(x2,a.x2); return *this;}
  const inline double8& operator-=(const double8 a) {x1 = _mm256_sub_pd(x1,a.x1); x2 = _mm256_sub_pd(x2,a.x2); return *this;}
  const inline double8& operator*=(const double8 a) {x1 = _mm256_mul_pd(x1,a.x1); x2 = _mm256_mul_pd(x2,a.x2); return *this;}
  const inline double8& operator/=(const double8 a) {x1 = _mm256_div_pd(x1,a.x1); x2 = _mm256_div_pd(x2,a.x2); return *this;}
  const inline double8 operator -()  const {return double8(_mm256_sub_pd(_mm256_set1_pd(0.0),x1),_mm256_sub_pd(_mm256_set1_pd(0.0),x2)); }

  friend inline const double8 operator + (double8 a, const double8 b) { return a+=b; };
  friend inline const double8 operator - (double8 a, const double8 b) { return a-=b; };
  friend inline const double8 operator * (double8 a, const double8 b) { return a*=b; };
  friend inline const double8 operator / (double8 a, const double8 b) { return a/=b; };

  /* basic comparison */
  const inline double8 operator!=(const double8 a) const { return double8(_mm256_cmp_pd(x1,a.x1,_CMP_NEQ_OQ),_mm256_cmp_pd(x2,a.x2,_CMP_NEQ_OQ)); }
  const inline double8 operator==(const double8 a) const { return double8(_mm256_cmp_pd(x1,a.x1,_CMP_EQ_OQ),_mm256_cmp_pd(x2,a.x2,_CMP_EQ_OQ)); }
  const inline double8 operator< (const double8 a) const { return double8(_mm256_cmp_pd(x1,a.x1,_CMP_LT_OQ),_mm256_cmp_pd(x2,a.x2,_CMP_LT_OQ)); }
  const inline double8 operator<=(const double8 a) const { return double8(_mm256_cmp_pd(x1,a.x1,_CMP_LE_OQ),_mm256_cmp_pd(x2,a.x2,_CMP_LE_OQ)); }
  const inline double8 operator> (const double8 a) const { return double8(_mm256_cmp_pd(x1,a.x1,_CMP_GT_OQ),_mm256_cmp_pd(x2,a.x2,_CMP_GT_OQ)); }
  const inline double8 operator>=(const double8 a) const { return double8(_mm256_cmp_pd(x1,a.x1,_CMP_GE_OQ),_mm256_cmp_pd(x2,a.x2,_CMP_GE_OQ)); }

  const inline double8 operator& (const double8 rhs) const { return double8(_mm256_and_pd(x1, rhs.x1), _mm256_and_pd(x2,rhs.x2)); } 
  const inline double8 operator&&(const double8 rhs) const { return double8(_mm256_and_pd(x1, rhs.x1), _mm256_and_pd(x2,rhs.x2)); } 
  const inline double8 operator| (const double8 rhs) const { return double8(_mm256_or_pd(x1, rhs.x1), _mm256_or_pd(x2,rhs.x2)); } 
  const inline double8 operator||(const double8 rhs) const { return double8(_mm256_or_pd(x1, rhs.x1), _mm256_or_pd(x2,rhs.x2)); } 
  const inline double8 operator^ (const double8 rhs) const { return double8(_mm256_xor_pd(x1, rhs.x1), _mm256_xor_pd(x2,rhs.x2)); } 
  const inline double8 operator! () const { return double8(_mm256_xor_pd(x1, full_mask()), _mm256_xor_pd(x2,full_mask())); } 

  friend inline const double8 operator& (const bool lhs, double8 rhs) { return double8(_mm256_and_pd(lhs ? full_mask() : empty_mask(), rhs.x1),_mm256_and_pd(lhs ? full_mask() : empty_mask(), rhs.x2)); } 
  friend inline const double8 operator&&(const bool lhs, double8 rhs) { return double8(_mm256_and_pd(lhs ? full_mask() : empty_mask(), rhs.x1),_mm256_and_pd(lhs ? full_mask() : empty_mask(), rhs.x2)); } 
  friend inline const double8 operator| (const bool lhs, double8 rhs) { return double8(_mm256_or_pd(lhs ? full_mask() : empty_mask(), rhs.x1),_mm256_or_pd(lhs ? full_mask() : empty_mask(), rhs.x2)); } 
  friend inline const double8 operator||(const bool lhs, double8 rhs) { return double8(_mm256_or_pd(lhs ? full_mask() : empty_mask(), rhs.x1),_mm256_or_pd(lhs ? full_mask() : empty_mask(), rhs.x2)); } 
  friend inline const double8 operator^ (const bool lhs, double8 rhs) { return double8(_mm256_xor_pd(lhs ? full_mask() : empty_mask(), rhs.x1),_mm256_xor_pd(lhs ? full_mask() : empty_mask(), rhs.x2)); } 

#if 0
  friend inline const double8 operator& (const bool lhs, double8 rhs) { return double8(lhs & rhs.x1, lhs & rhs.x2); }
  friend inline const double8 operator&&(const bool lhs, double8 rhs) { return double8(lhs && rhs.x1, lhs && rhs.x2); }
  friend inline const double8 operator| (const bool lhs, double8 rhs) { return double8(lhs | rhs.x1, lhs | rhs.x2); }
  friend inline const double8 operator||(const bool lhs, double8 rhs) { return double8(lhs || rhs.x1, lhs || rhs.x2); }
  friend inline const double8 operator^ (const bool lhs, double8 rhs) { return double8(lhs ^ rhs.x1, lhs ^ rhs.x2); }
#endif

};


/***** LOAD/STORE *******/


/* vector load */
template<> static inline double8 vload<double8>(const double *data) 
{ return double8(data); }
/* vector store */
template<> static inline double8 vstore<double8>(double *ptr, const double8 data) 
{  _mm256_storeu_pd(ptr,data.x1), _mm256_storeu_pd(ptr+4,data.x2); }
/* gather loads */
template<> static inline double8 vgather<double8,int>(const double *data, const int addr) 
{  return double8(data[addr]);  }
template<> static inline double8 vgather<double8,int8>(const double *data, const int8 addr) 
{  return double8(data,addr);  }
/* scatter store */
template<> static inline double8 vscatter<double8,int8>(double *data, const int8 addr, const double8 value)
{ 
  assert(0);
  return value;
}

/******** MATH ***********/

static inline double8 __min    (const double8 a, const double8 b) {return double8(_mm256_min_pd(a.x1,b.x1), _mm256_min_pd(a.x2,b.x2)); }
static inline double8 __max    (const double8 a, const double8 b) {return double8(_mm256_max_pd(a.x1,b.x1), _mm256_max_pd(a.x2,b.x2)); }
static inline double8 __abs    (const double8 x)                  {return double8(_mm256_and_pd(x.x1, abs_mask()),_mm256_and_pd(x.x2,abs_mask()));  }
static inline double8 __pow    (const double8 x, const double8 p) {return double8(_mm256_pow_pd(x.x1,p.x1),_mm256_pow_pd(x.x2,p.x2));}
static inline double8 __exp    (const double8 x)                  {return double8(_mm256_exp_pd(x.x1),_mm256_exp_pd(x.x2));}
static inline double8 __log    (const double8 x)                  {return double8(_mm256_log_pd(x.x1),_mm256_log_pd(x.x2));}
static inline double8 __log10  (const double8 x)                  {return double8(0.434294481903252) * __log(x);}
static inline double8 __sqrt   (const double8 x)                  {return double8(_mm256_sqrt_pd(x.x1),_mm256_sqrt_pd(x.x2)); }
static inline double8 __invsqrt(const double8 x)                  {return double8(_mm256_invsqrt_pd(x.x1),_mm256_invsqrt_pd(x.x2)); }
static inline double8 __cbrt   (const double8 x)                  {return double8(_mm256_cbrt_pd(x.x1),_mm256_cbrt_pd(x.x2)); }
static inline double8 __invcbrt(const double8 x)                  {return double8(_mm256_invcbrt_pd(x.x1),_mm256_invcbrt_pd(x.x2));  }
static inline double8 __cos    (const double8 x)                  {return double8(_mm256_cos_pd(x.x1),_mm256_cos_pd(x.x2));}
static inline double8 __sin    (const double8 x)                  {return double8(_mm256_sin_pd(x.x1),_mm256_sin_pd(x.x2));}

/**********************************/
/********* CONDITIONALS ***********/
/**********************************/

template<> static inline double8 vtrue <double8>() { return double8(full_mask(), full_mask()); }
template<> static inline double8 vfalse <double8>() { return double8(empty_mask(), empty_mask()); }

template<> struct Mask<double8>
{
  int imask;
  double8 mask;
  explicit inline Mask(const double8 condition) : 
    imask(_mm256_movemask_pd(condition.x1) | (_mm256_movemask_pd(condition.x2)<<4)), mask(condition) {}
  inline operator double8() const { return mask; }
  explicit inline Mask(const bool condition) : imask(condition ? 0xFF : 0x0), mask(condition ? vtrue<double8>() : vfalse<double8>()) {}
  inline void operator() (double8 &a, const double8 b) const __attribute__((always_inline))
  {
    a.x1 = _mm256_blendv_pd(a.x1,b.x1,mask.x1);
    a.x2 = _mm256_blendv_pd(a.x2,b.x2,mask.x2);
  }
  inline bool any() const { return imask; }
  inline bool all() const { return imask == 0xFF; }
  inline bool none() const {return imask == 0; }
  inline void inverse() { imask ^= 0xFF; mask = double8(mask) ^ double8(vtrue<double8>());}
  inline Mask operator&&(const Mask m) const {return Mask(m.mask && mask); }
};


typedef double8 vdouble;
typedef int8    vinteger;
typedef double8 vboolean;
