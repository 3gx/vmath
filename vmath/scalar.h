#pragma once
#include <cmath>
#include <algorithm>
#include "vmath.h"

/* scalar implementation of vmath */

/************************/
/****** INTEGER *********/
/************************/

template<> static inline int vspan<int>(const int i) { return i; }

/***** LOAD/STORE *******/

/******** MATH ***********/

static inline int __max(const int a, const int b) { return std::max(a,b); }
static inline int __min(const int a, const int b) { return std::min(a,b); }
static inline int __abs(const int a)              { return std::abs(a); }

/************************/
/******* DOUBLE *********/
/************************/

/***** LOAD/STORE *******/

/* vector load */
template<> static inline double  vload<double >(const double *data) 
{ return *data; } 
/* vector store */
template<> static inline double vstore<double>(double *ptr, const double data) 
{   *ptr = data;  return data; }
/* gather load */
template<> static inline double vgather<double,int>(const double *data, const int addr) 
{  return data[addr]; }
/* scatter stores */
template<> static inline double vscatter<double,int>(double *data, int addr, const double value) 
{ data[addr] = value; return value; }

/******** MATH ***********/

static inline double __max    (const double a, const double b) { return std::max(a,b); }
static inline double __min    (const double a, const double b) { return std::min(a,b); }
static inline double __abs    (const double x)                 { return std::abs(x); }
static inline double __pow    (const double x, const double p) { return std::pow(x,p); }
static inline double __exp    (const double x)                 { return std::exp(x); }
static inline double __log    (const double x)                 { return std::log(x); }
static inline double __log10  (const double x)                 { return 0.434294481903252 * __log(x);}
static inline double __sqrt   (const double x)                 { return std::sqrt(x); }
static inline double __invsqrt(const double x)                 { return 1.0/std::sqrt(x); }
static inline double __cbrt   (const double x)                 { return std::pow(x, 1.0/3.0); }
static inline double __invcbrt(const double x)                 { return std::pow(x, -1.0/3.0); }
static inline double __cos    (const double x)                 { return std::cos(x); }
static inline double __sin    (const double x)                 { return std::sin(x); }

/**********************************/
/********* CONDITIONALS ***********/
/**********************************/

template<> static inline bool vtrue <bool>() { return true;  }
template<> static inline bool vfalse<bool>() { return false; }
template<> struct Mask<bool>
{
  bool mask;
  explicit inline Mask(const  bool  condition) : mask(condition      ) {}
  inline operator bool() const { return mask; }
  template<typename T>
    inline void operator() (T &a, const T b) const __attribute__((always_inline))
    {
      if (mask) a = b;
    }
  inline bool any() const {return mask;}
  inline bool all() const {return mask;}
  inline bool none() const {return !mask;}
  inline void inverse() { mask = !mask; }
  inline Mask operator&&(const Mask _mask) const {return Mask(_mask.mask && mask); }
};

