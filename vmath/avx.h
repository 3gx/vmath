#include <immintrin.h>
#include <iostream>

#include "vmath.h"

/* avx implementation of vmath */

#ifndef __AVX__
#error "AVX support must be enabled!"
#endif

/************************/
/****** INTEGER *********/
/************************/

/* vector integer */
struct int4
{
  enum {SIZE = 4};
  __m128i i;

  inline int4() {}
  inline int4(const __m128i v) : i(v) {}
  inline int4(const int  s) : i(_mm_set1_epi32(s)) {}
  inline int4(const int *data, const int idx[4]) : i(_mm_set_epi32(data[idx[3]], data[idx[2]], data[idx[1]], data[idx[0]])) {}
  inline operator __m128i() const {return i;}
  inline operator __m256d() const {return _mm256_cvtepi32_pd(i);}

  const inline int4& operator =(const int4 a) {i  = a.i; return *this;}
  const inline int4& operator =(const __m128i a) {i  = a; return *this;}

  const inline int4& operator+=(const int4 a) {i = _mm_add_epi32(i,a.i); return *this;}
  const inline int4& operator-=(const int4 a) {i = _mm_sub_epi32(i,a.i); return *this;}
  const inline int4& operator*=(const int4 a) {i = _mm_mullo_epi32(i,a.i); return *this;}
  const inline int4& operator/=(const int4 a) {i = _mm_div_epi32(i,a.i); return *this;}

  friend inline int4 operator + (int4 a, const int4 b) { return a+=b; };
  friend inline int4 operator - (int4 a, const int4 b) { return a-=b; };
  friend inline int4 operator * (int4 a, const int4 b) { return a*=b; };
  friend inline int4 operator / (int4 a, const int4 b) { return a/=b; };

  inline int operator[] (const int idx) const
  {
    union {__m128i v; int s[4];} val = {i};
    return val.s[idx];
  }

  inline friend std::ostream &operator << (std::ostream &ofs, const int4 v)
  {
    ofs << "[ " <<
      v[0] << " " << 
      v[1] << " " << 
      v[2] << " " << 
      v[3] << " ]";
    return ofs;
  }
};

template<> static inline int4 vspan<int4>(const int i) { return _mm_set_epi32(i+3,i+2,i+1,i); }

/***** LOAD/STORE *******/

/******** MATH ***********/

static inline int4 __max(const int4 a, const int4 b) { return int4(_mm_max_epi32(a,b)); }
static inline int4 __min(const int4 a, const int4 b) { return int4(_mm_min_epi32(a,b)); }
static inline int4 __abs(const int4 a)               { return int4(_mm_abs_epi32(a)); }

/************************/
/******* DOUBLE *********/
/************************/

static inline __m256d get_mask(const unsigned long long mask) { return(__m256d)_mm256_set1_epi64x(mask); }
static inline __m256d abs_mask() { return get_mask(~(1ULL<<63)); }
static inline __m256d full_mask() { return get_mask(-1); }
static inline __m256d empty_mask() { return get_mask(0); }

/* vector double */
struct double4
{
  enum {SIZE = 4};
  inline static int round(const int n) { return n & (-SIZE); }


  __m256d x;
  inline double4() {}
  inline double4(const __m256d v) : x(v) {}
  inline double4(const double  s) : x(_mm256_set1_pd(s)) {}
  inline double4(const double *data, const int idx[4])
  {
    x = _mm256_set_pd(data[idx[3]], data[idx[2]], data[idx[1]], data[idx[0]]);
  }
  inline double4(const double *data, const int4 _idx)
  {
    union { __m128i v; int s[4]; } idx = {_idx};
    x = double4(data, idx.s);
  }
  inline double4(const double *data) : x(_mm256_loadu_pd(data)) {}

  inline operator __m256d() const {return x;}
  inline operator __m128i() const {return _mm256_cvttpd_epi32(x);}

  /* basic arithmetics */
  const inline double4& operator =(const double4 a) {x  = a.x; return *this;}
  const inline double4& operator =(const __m256d a) {x  = a; return *this;}
  const inline double4& operator =(const double  a) {x = _mm256_set1_pd(a); return *this;}
  const inline double4& operator+=(const double4 a) {x = _mm256_add_pd(x,a.x); return *this;}
  const inline double4& operator-=(const double4 a) {x = _mm256_sub_pd(x,a.x); return *this;}
  const inline double4& operator*=(const double4 a) {x = _mm256_mul_pd(x,a.x); return *this;}
  const inline double4& operator/=(const double4 a) {x = _mm256_div_pd(x,a.x); return *this;}
  const inline double4 operator -()  const {return double4(_mm256_sub_pd(double4(0.0), x));}

  friend inline const double4 operator + (double4 a, const double4 b) { return a+=b; };
  friend inline const double4 operator - (double4 a, const double4 b) { return a-=b; };
  friend inline const double4 operator * (double4 a, const double4 b) { return a*=b; };
  friend inline const double4 operator / (double4 a, const double4 b) { return a/=b; };

  /* basic comparison */
  const inline double4 operator!=(const double4 a) const { return double4(_mm256_cmp_pd(x,a,_CMP_NEQ_OQ)); }
  const inline double4 operator==(const double4 a) const { return double4(_mm256_cmp_pd(x,a,_CMP_EQ_OQ)); }
  const inline double4 operator< (const double4 a) const { return double4(_mm256_cmp_pd(x,a,_CMP_LT_OQ)); }
  const inline double4 operator<=(const double4 a) const { return double4(_mm256_cmp_pd(x,a,_CMP_LE_OQ)); }
  const inline double4 operator> (const double4 a) const { return double4(_mm256_cmp_pd(x,a,_CMP_GT_OQ)); }
  const inline double4 operator>=(const double4 a) const { return double4(_mm256_cmp_pd(x,a,_CMP_GE_OQ)); }

  const inline double4 operator& (const double4 rhs) const { return double4(_mm256_and_pd(x, rhs)); } 
  const inline double4 operator&&(const double4 rhs) const { return double4(_mm256_and_pd(x, rhs)); } 
  const inline double4 operator| (const double4 rhs) const { return double4(_mm256_or_pd (x, rhs)); } 
  const inline double4 operator||(const double4 rhs) const { return double4(_mm256_or_pd (x, rhs)); } 
  const inline double4 operator^ (const double4 rhs) const { return double4(_mm256_xor_pd (x, rhs)); } 
  const inline double4 operator! () const { return double4(_mm256_xor_pd(x, full_mask())); }

  friend inline const double4 operator& (const bool lhs, double4 rhs) { return double4(_mm256_and_pd(lhs ? full_mask() : empty_mask(), rhs)); } 
  friend inline const double4 operator&&(const bool lhs, double4 rhs) { return double4(_mm256_and_pd(lhs ? full_mask() : empty_mask(), rhs)); }
  friend inline const double4 operator| (const bool lhs, double4 rhs) { return double4(_mm256_or_pd (lhs ? full_mask() : empty_mask(), rhs)); }
  friend inline const double4 operator||(const bool lhs, double4 rhs) { return double4(_mm256_or_pd (lhs ? full_mask() : empty_mask(), rhs)); }
  friend inline const double4 operator^ (const bool lhs, double4 rhs) { return double4(_mm256_xor_pd(lhs ? full_mask() : empty_mask(), rhs)); }


  inline double operator[] (const int i) const
  {
    union {__m256d v; double s[4];} val = {x};
    return val.s[i];
  }
  inline friend std::ostream &operator << (std::ostream &ofs, const double4 v)
  {
    ofs << "[ " <<
      v[0] << " " << 
      v[1] << " " << 
      v[2] << " " << 
      v[3] << " ]";
    return ofs;
  }
};


/***** LOAD/STORE *******/


/* vector load */
template<> static inline double4 vload<double4>(const double *data) 
{ return double4(_mm256_loadu_pd(data)); }
/* vector store */
template<> static inline double4 vstore<double4>(double *ptr, const double4 data) 
{  _mm256_storeu_pd(ptr,data);   return data; }
/* gather loads */
template<> static inline double4 vgather<double4,int>(const double *data, const int addr) 
{  return double4(data[addr]);  }
template<> static inline double4 vgather<double4,int4>(const double *data, const int4 addr) 
{  return double4(data,addr);  }
/* scatter store */
template<> static inline double4 vscatter<double4,int4>(double *data, const int4 addr, const double4 value)
{ 
  union {__m128i v; int    s[4];} idx = {addr};
  union {__m256d v; double s[4];} val = {value};
  data[idx.s[3]] = val.s[3];
  data[idx.s[2]] = val.s[2];
  data[idx.s[1]] = val.s[1];
  data[idx.s[0]] = val.s[0];
  return value;
}

/******** MATH ***********/

static inline double4 __min    (const double4 a, const double4 b) {return double4(_mm256_min_pd(a,b)); }
static inline double4 __max    (const double4 a, const double4 b) {return double4(_mm256_max_pd(a,b)); }
static inline double4 __abs    (const double4 x)                  {return double4(_mm256_and_pd(x, abs_mask()));  }
static inline double4 __pow    (const double4 x, const double4 p) {return double4(_mm256_pow_pd(x,p));}
static inline double4 __exp    (const double4 x)                  {return double4(_mm256_exp_pd(x));}
static inline double4 __log    (const double4 x)                  {return double4(_mm256_log_pd(x));}
static inline double4 __log10  (const double4 x)                  {return double4(0.434294481903252) * __log(x);}
static inline double4 __sqrt   (const double4 x)                  {return double4(_mm256_sqrt_pd(x)); }
static inline double4 __invsqrt(const double4 x)                  {return double4(_mm256_invsqrt_pd(x)); }
static inline double4 __cbrt   (const double4 x)                  {return double4(_mm256_cbrt_pd(x)); }
static inline double4 __invcbrt(const double4 x)                  {return double4(_mm256_invcbrt_pd(x));  }
static inline double4 __cos    (const double4 x)                  {return double4(_mm256_cos_pd(x));}
static inline double4 __sin    (const double4 x)                  {return double4(_mm256_sin_pd(x));}

/**********************************/
/********* CONDITIONALS ***********/
/**********************************/

template<> static inline double4 vtrue <double4>() { return double4(full_mask());  }
template<> static inline double4 vfalse<double4>() { return double4(empty_mask()); }

template<> struct Mask<double4>
{
  int imask;
  double4 mask;
  explicit inline Mask(const double4 condition) : imask(_mm256_movemask_pd(condition)), mask(condition) {}
  explicit inline Mask(const bool condition) : imask(condition ? 0xF : 0x0), mask(condition ? full_mask() : empty_mask()) {}
  inline operator __m256d() const { return mask; }
  inline void operator() (double4 &a, const double4 b) const __attribute__((always_inline))
  {
    a = _mm256_blendv_pd(a,b,mask);
  }
  inline bool any() const { return imask; }
  inline bool all() const { return imask == 0xF; }
  inline bool none() const {return imask == 0; }
  inline void inverse() { imask ^= 0xF; mask = double4(mask) ^ double4(full_mask());}
  inline Mask operator&&(const Mask m) const {return Mask(m.mask && mask); }
};


typedef double4 vdouble;
typedef int4    vinteger;
typedef double4 vboolean;
