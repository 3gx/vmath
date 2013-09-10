#include <immintrin.h>
#include <iostream>

#include "vmath.h"

/* avx implementation of vmath */

#ifndef __SSE__
#error "SSE support must be enabled!"
#endif

/************************/
/****** INTEGER *********/
/************************/

/* vector integer */
struct int2
{
  enum {SIZE = 2};
  __m128i i;

  inline int2() {}
  inline int2(const __m128i v) : i(v) {}
  inline int2(const int  s) : i(_mm_set1_epi32(s)) {}
  inline int2(const int *data, const int idx[2]) : i(_mm_set_epi32(data[idx[1]], data[idx[1]], data[idx[1]], data[idx[0]])) {}
  inline operator __m128i() const {return i;}
  inline operator __m128d() const {return _mm_cvtepi32_pd(i);}

  const inline int2& operator =(const int2 a) {i  = a.i; return *this;}
  const inline int2& operator =(const __m128i a) {i  = a; return *this;}

  const inline int2& operator+=(const int2 a) {i = _mm_add_epi32(i,a.i); return *this;}
  const inline int2& operator-=(const int2 a) {i = _mm_sub_epi32(i,a.i); return *this;}
  const inline int2& operator*=(const int2 a) {i = _mm_mullo_epi32(i,a.i); return *this;}
  const inline int2& operator/=(const int2 a) {i = _mm_div_epi32(i,a.i); return *this;}

  friend inline int2 operator + (int2 a, const int2 b) { return a+=b; };
  friend inline int2 operator - (int2 a, const int2 b) { return a-=b; };
  friend inline int2 operator * (int2 a, const int2 b) { return a*=b; };
  friend inline int2 operator / (int2 a, const int2 b) { return a/=b; };

  inline int operator[] (const int idx) const
  {
    union {__m128i v; int s[4];} val = {i};
    return val.s[idx];
  }

  inline friend std::ostream &operator << (std::ostream &ofs, const int2 v)
  {
    ofs << "[ " <<
      v[0] << " " << 
      v[1] << " ]";
    return ofs;
  }
};

/***** LOAD/STORE *******/

/******** MATH ***********/

static inline int2 __max(const int2 a, const int2 b) { return int2(_mm_max_epi32(a,b)); }
static inline int2 __min(const int2 a, const int2 b) { return int2(_mm_min_epi32(a,b)); }
static inline int2 __abs(const int2 a)               { return int2(_mm_abs_epi32(a)); }

/************************/
/******* DOUBLE *********/
/************************/

static inline __m128d get_mask(const unsigned long long mask) { return(__m128d)_mm_set1_epi64x(mask); }
static inline __m128d abs_mask() { return get_mask(~(1ULL<<63)); }
static inline __m128d full_mask() { return get_mask(-1); }
static inline __m128d empty_mask() { return get_mask(0); }

/* vector double */
struct double2
{
  enum {SIZE = 2};
  inline static int round(const int n) { return n & (-SIZE); }


  __m128d x;
  inline double2() {}
  inline double2(const __m128d v) : x(v) {}
  inline double2(const double  s) : x(_mm_set1_pd(s)) {}
  inline double2(const double *data, const int idx[2])
  {
    x = _mm_set_pd(data[idx[1]], data[idx[0]]);
  }
  inline double2(const double *data, const int2 _idx)
  {
    union { __m128i v; int s[4]; } idx = {_idx};
    const int __idx[2] = {idx.s[0], idx.s[1]};
    x = double2(data, __idx);
  }
  inline double2(const double *data) : x(_mm_loadu_pd(data)) {}

  inline operator __m128d() const {return x;}
  inline operator __m128i() const {return _mm_cvttpd_epi32(x);}

  /* basic arithmetics */
  const inline double2& operator =(const double2 a) {x  = a.x; return *this;}
  const inline double2& operator =(const __m128d a) {x  = a; return *this;}
  const inline double2& operator =(const double  a) {x = _mm_set1_pd(a); return *this;}
  const inline double2& operator+=(const double2 a) {x = _mm_add_pd(x,a.x); return *this;}
  const inline double2& operator-=(const double2 a) {x = _mm_sub_pd(x,a.x); return *this;}
  const inline double2& operator*=(const double2 a) {x = _mm_mul_pd(x,a.x); return *this;}
  const inline double2& operator/=(const double2 a) {x = _mm_div_pd(x,a.x); return *this;}
  const inline double2 operator -()  const {return double2(_mm_sub_pd(double2(0.0), x));}

  friend inline const double2 operator + (double2 a, const double2 b) { return a+=b; };
  friend inline const double2 operator - (double2 a, const double2 b) { return a-=b; };
  friend inline const double2 operator * (double2 a, const double2 b) { return a*=b; };
  friend inline const double2 operator / (double2 a, const double2 b) { return a/=b; };

  /* basic comparison */
  const inline double2 operator!=(const double2 a) const { return double2(_mm_cmpneq_pd(x,a)); }
  const inline double2 operator==(const double2 a) const { return double2(_mm_cmpeq_pd (x,a)); }
  const inline double2 operator< (const double2 a) const { return double2(_mm_cmplt_pd (x,a)); }
  const inline double2 operator<=(const double2 a) const { return double2(_mm_cmple_pd (x,a)); }
  const inline double2 operator> (const double2 a) const { return double2(_mm_cmpgt_pd (x,a)); }
  const inline double2 operator>=(const double2 a) const { return double2(_mm_cmpge_pd (x,a)); }

  const inline double2 operator& (const double2 rhs) const { return double2(_mm_and_pd(x, rhs)); } 
  const inline double2 operator&&(const double2 rhs) const { return double2(_mm_and_pd(x, rhs)); } 
  const inline double2 operator| (const double2 rhs) const { return double2(_mm_or_pd (x, rhs)); } 
  const inline double2 operator||(const double2 rhs) const { return double2(_mm_or_pd (x, rhs)); } 
  const inline double2 operator^ (const double2 rhs) const { return double2(_mm_xor_pd (x, rhs)); } 
  const inline double2 operator! () const { return double2(_mm_xor_pd(x, full_mask())); }

  friend inline const double2 operator& (const bool lhs, double2 rhs) { return double2(_mm_and_pd(lhs ? full_mask() : empty_mask(), rhs)); } 
  friend inline const double2 operator&&(const bool lhs, double2 rhs) { return double2(_mm_and_pd(lhs ? full_mask() : empty_mask(), rhs)); }
  friend inline const double2 operator| (const bool lhs, double2 rhs) { return double2(_mm_or_pd (lhs ? full_mask() : empty_mask(), rhs)); }
  friend inline const double2 operator||(const bool lhs, double2 rhs) { return double2(_mm_or_pd (lhs ? full_mask() : empty_mask(), rhs)); }
  friend inline const double2 operator^ (const bool lhs, double2 rhs) { return double2(_mm_xor_pd(lhs ? full_mask() : empty_mask(), rhs)); }


  inline double operator[] (const int i) const
  {
    union {__m128d v; double s[4];} val = {x};
    return val.s[i];
  }
  inline friend std::ostream &operator << (std::ostream &ofs, const double2 v)
  {
    ofs << "[ " <<
      v[0] << " " << 
      v[1] << " ]";
    return ofs;
  }
};


/***** LOAD/STORE *******/


/* vector load */
template<> static inline double2 vload<double2>(const double *data) 
{ return double2(_mm_loadu_pd(data)); }
/* vector store */
template<> static inline double2 vstore<double2>(double *ptr, const double2 data) 
{  _mm_storeu_pd(ptr,data);   return data; }
/* gather loads */
template<> static inline double2 vgather<double2,int>(const double *data, const int addr) 
{  return double2(data[addr]);  }
template<> static inline double2 vgather<double2,int2>(const double *data, const int2 addr) 
{  return double2(data,addr);  }
/* scatter store */
template<> static inline double2 vscatter<double2,int2>(double *data, const int2 addr, const double2 value)
{ 
  union {__m128i v; int    s[4];} idx = {addr};
  union {__m128d v; double s[2];} val = {value};
  data[idx.s[1]] = val.s[1];
  data[idx.s[0]] = val.s[0];
  return value;
}

/******** MATH ***********/

static inline double2 __min    (const double2 a, const double2 b) {return double2(_mm_min_pd(a,b)); }
static inline double2 __max    (const double2 a, const double2 b) {return double2(_mm_max_pd(a,b)); }
static inline double2 __abs    (const double2 x)                  {return double2(_mm_and_pd(x, abs_mask()));  }
static inline double2 __pow    (const double2 x, const double2 p) {return double2(_mm_pow_pd(x,p));}
static inline double2 __exp    (const double2 x)                  {return double2(_mm_exp_pd(x));}
static inline double2 __log    (const double2 x)                  {return double2(_mm_log_pd(x));}
static inline double2 __log10  (const double2 x)                  {return double2(0.434294481903252) * __log(x);}
static inline double2 __sqrt   (const double2 x)                  {return double2(_mm_sqrt_pd(x)); }
static inline double2 __invsqrt(const double2 x)                  {return double2(_mm_invsqrt_pd(x)); }
static inline double2 __cbrt   (const double2 x)                  {return double2(_mm_cbrt_pd(x)); }
static inline double2 __invcbrt(const double2 x)                  {return double2(_mm_invcbrt_pd(x));  }
static inline double2 __cos    (const double2 x)                  {return double2(_mm_cos_pd(x));}
static inline double2 __sin    (const double2 x)                  {return double2(_mm_sin_pd(x));}

/**********************************/
/********* CONDITIONALS ***********/
/**********************************/

template<> static inline double2 vtrue <double2>() { return double2(full_mask());  }
template<> static inline double2 vfalse<double2>() { return double2(empty_mask()); }

template<> struct Mask<double2>
{
  int imask;
  double2 mask;
  explicit inline Mask(const double2 condition) : imask(_mm_movemask_pd(condition)), mask(condition) {}
  explicit inline Mask(const bool condition) : imask(condition ? 0x3 : 0x0), mask(condition ? full_mask() : empty_mask()) {}
  inline operator __m128d() const { return mask; }
  inline void operator() (double2 &a, const double2 b) const __attribute__((always_inline))
  {
    a = _mm_blendv_pd(a,b,mask);
  }
  inline bool any() const { return imask; }
  inline bool all() const { return imask == 0x3; }
  inline bool none() const {return imask == 0; }
  inline void inverse() { imask ^= 0x3; mask = double2(mask) ^ double2(full_mask());}
  inline Mask operator&&(const Mask m) const {return Mask(m.mask && mask); }
};


typedef double2 vdouble;
typedef int2    vinteger;
typedef double2 vboolean;
