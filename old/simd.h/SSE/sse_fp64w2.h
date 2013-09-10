#pragma once

#ifndef __SSE__
#error "Please enable SSE support"
#endif

#define __out

#include <cassert>

namespace SIMD
{
  enum {WIDTH = 4};
  typedef double _v2df  __attribute__((vector_size(16)));
  typedef float  _v4sf  __attribute__((vector_size(16)));

  struct dvec 
  {
    _v2df v[2];
    dvec() {}
    public:
    dvec(_v2df a, const _v2df b) : v{a,b} {}
    dvec(_v2df a) : v{a,a} {}
    dvec(const double  f) : v{(_v2df){f,f}, (_v2df){f,f}} {}
    dvec(const double  a, const double b, const double c, const double d) : 
      v{(_v2df){a,b},(_v2df){c,d}} {}
    dvec(const double *p, const bool aligned = true) 
    {
      if (aligned)
      {
          v[0] = *(_v2df*)p;
          v[1] = *(_v2df*)(p+WIDTH);
      }
      else
      {
         v[0] =  __builtin_ia32_loadupd(p      );
         v[1] =  __builtin_ia32_loadupd(p+WIDTH);
      }
    }
    
    _v2df operator[](const int i) const {return v[i];}

		const dvec& operator =(const dvec a) {v[0]  = a.v[0]; v[1]  = a.v[1]; return *this;}
		const dvec& operator+=(const dvec a) {v[0] += a.v[0]; v[1] += a.v[1]; return *this;}
		const dvec& operator-=(const dvec a) {v[0] -= a.v[0]; v[1] -= a.v[1]; return *this;}
		const dvec& operator*=(const dvec a) {v[0] *= a.v[0]; v[1] *= a.v[1]; return *this;}
		const dvec& operator/=(const dvec a) {v[0] /= a.v[0]; v[1] /= a.v[1]; return *this;}
		const dvec operator -()  const {return dvec(-v[0], -v[1]);}
		const dvec operator +(const dvec a) const {return dvec(v[0] + a.v[0], v[1] + a.v[1]);}
		const dvec operator -(const dvec a) const {return dvec(v[0] - a.v[0], v[1] - a.v[1]);}
		const dvec operator *(const dvec a) const {return dvec(v[0] * a.v[0], v[1] * a.v[1]);}
		const dvec operator /(const dvec a) const {return dvec(v[0] / a.v[0], v[1] / a.v[1]);}
		
    const dvec operator!=(const dvec &rhs) const 
    {
			return dvec(
          __builtin_ia32_cmpneqpd(v[0], rhs.v[0]),
          __builtin_ia32_cmpneqpd(v[1], rhs.v[1]));
		}
		const dvec operator==(const dvec rhs) const {
			return dvec(
          __builtin_ia32_cmpeqpd(v[0], rhs.v[0]),
          __builtin_ia32_cmpeqpd(v[1], rhs.v[1]));
		}

		const dvec operator<(const dvec &rhs) const{
			return dvec(
          __builtin_ia32_cmpltpd(v[0], rhs.v[0]),
          __builtin_ia32_cmpltpd(v[1], rhs.v[1]));
		}
		const dvec operator<=(const dvec &rhs) const{
			return dvec(
          __builtin_ia32_cmplepd(v[0], rhs.v[0]),
          __builtin_ia32_cmplepd(v[1], rhs.v[1]));
		}
		const dvec operator>(const dvec &rhs) const{
			return dvec(
          __builtin_ia32_cmpgtpd(v[0], rhs.v[0]),
          __builtin_ia32_cmpgtpd(v[1], rhs.v[1]));
		}
		const dvec operator>=(const dvec &rhs) const{
			return dvec(
          __builtin_ia32_cmpgepd(v[0], rhs.v[0]),
          __builtin_ia32_cmpgepd(v[1], rhs.v[1]));
		}
		const dvec operator|(const dvec &rhs) const{
			return dvec(
          __builtin_ia32_orpd(v[0], rhs.v[0]),
          __builtin_ia32_orpd(v[1], rhs.v[1]));
		}
		const dvec operator&(const dvec &rhs) const{
			return dvec(
          __builtin_ia32_andpd(v[0], rhs.v[0]),
          __builtin_ia32_andpd(v[1], rhs.v[1]));
		}

    operator bool() const {return true;}
		
    const double ch(const int i) const {
			switch(i) {
				case 0:	return __builtin_ia32_vec_ext_v2df(v[0], 0);
				case 1:	return __builtin_ia32_vec_ext_v2df(v[0], 1);
				case 2:	return __builtin_ia32_vec_ext_v2df(v[1], 0);
				case 3:	return __builtin_ia32_vec_ext_v2df(v[1], 1);
				default: assert(0);
			}
    }
	

  };

#if 0
  inline static void transpose(const dvec v0, const dvec v1, dvec &t0, dvec &t1)
  {
    t0 = __builtin_ia32_unpcklpd(v0, v1);
    t1 = __builtin_ia32_unpckhpd(v0, v1);
  }

  template<const int N>
    inline static void transpose(dvec *d)
    {
      dvec v1, v2;
      switch(N)
      {
        case 1:
          transpose(d[0], d[1], d[0], d[1]);
          break;
        case 2:
          transpose(d[0], d[2], d[0], d[2]);
          transpose(d[5], d[7], d[5], d[7]);
          transpose(d[1], d[3], v1, v2);
          transpose(d[4], d[6], d[1], d[3]);
          d[4] = v1;
          d[6] = v2;
          break;
        default:
          assert(N == 1);
      }
    }

  /* slow */
  template<const int NX, const int NY>
    static void transpose(const double* __restrict__ _src, double* __restrict__ _dst)
    {
      assert((NX%SIMD::WIDTH) == 0);
      assert((NY%SIMD::WIDTH) == 0);

      const int M = NX/SIMD::WIDTH;
      const int N = NY/SIMD::WIDTH;
      const dvec (*src)[M] = (const dvec (*)[M])_src;
      __out dvec (*dst)[M] = (      dvec (*)[M])_dst;

      for (int j = 0; j < N; j++)
      {
        const int j2 = j << 1;
        __builtin_prefetch(&src[j2+2][0]);
        __builtin_prefetch(&src[j2+3][0]);
        for (int i = 0; i < M; i++)
        {
          int i2 = i << 1;
          __builtin_prefetch(&dst[i2+2][j]);
          __builtin_prefetch(&dst[i2+3][j]);
          __builtin_prefetch(&dst[i2+4][j]);
          __builtin_prefetch(&dst[i2+5][j]);
          __builtin_prefetch(&dst[i2+6][j]);
          __builtin_prefetch(&dst[i2+7][j]);
          dvec a = src[j2  ][i];
          dvec b = src[j2+1][i];
          transpose(a, b, a, b);
          dst[i2  ][j] = a;
          dst[i2+1][j] = b;
        }
      }

#if 0
      for (int j = 0; j < NY; j++)
        for (int i = 0; i < NX; i++)
          assert(_src[j*NX+i] == _dst[i*NX+j]);
#endif
    }

  template<class T>
    inline static void sstore(T &ptr, const dvec val)
    {
      __builtin_ia32_movntpd((double*)&ptr, val);
    }

  template<const int N>
    inline static dvec broadcast(const dvec x)
    {
      const int mask =
        N == 0 ? (0 << 0) + (0 << 0) : (1 << 0) + (1 << 1);
      return __builtin_ia32_shufpd(x, x, mask);
    }
#endif
  
  template<const bool ALIGN, typename T>
    inline static dvec store(T &ptr, const dvec val)
    {
      if (ALIGN)
      {
        _v2df* vptr = (_v2df*)&ptr;
        vptr[0] = val[0];
        vptr[1] = val[1];
      }
      else
      {
        __builtin_ia32_storeupd( (double*)&ptr,          val[0]);
        __builtin_ia32_storeupd(((double*)&ptr) + WIDTH, val[1]);
      }
      return val;
    }
  
  template<const bool ALIGN, typename T>
    inline static dvec load(const T &ptr)
    {
      return dvec((double*)&ptr, ALIGN);
    }

  
  inline static dvec sqrt(const dvec x)
  {
    return dvec(
        __builtin_ia32_sqrtpd(x[0]),
        __builtin_ia32_sqrtpd(x[1]));
  }
  
  inline static dvec rsqrt(const dvec x)
  {
#if 0
    const fvec y = __builtin_ia32_rsqrtps(x);
    return (fvec(-0.5f) * y) * (x*y*y + fvec(-3.0f));
#else
    const _v4sf y = __builtin_ia32_movlhps(__builtin_ia32_cvtpd2ps(x[0]),__builtin_ia32_cvtpd2ps(x[1]));
    const _v4sf ysp = __builtin_ia32_rsqrtps(y);
    const _v2df y1 = __builtin_ia32_cvtps2pd(ysp);
    const _v2df y2 = __builtin_ia32_cvtps2pd(__builtin_ia32_movhlps(ysp, ysp));
    return dvec(y1, y2);
#if 0
    const v2df c1 = {-0.5, -0.5};
    const v2df c2 = {-3.0, -3.0};
    y1 = (c1 * y1) * (r2_1*y1*y1 + c2);
    y2 = (c1 * y2) * (r2_2*y2*y2 + c2);
    y1 = (c1 * y1) * (r2_1*y1*y1 + c2);
    y2 = (c1 * y2) * (r2_2*y2*y2 + c2);
    const v2df flag_1  = __builtin_ia32_cmpgtpd(r2_1, (v2df){0.0, 0.0});
    const v2df flag_2  = __builtin_ia32_cmpgtpd(r2_2, (v2df){0.0, 0.0});
    r2_1  = __builtin_ia32_andpd(flag_1, y1);
    r2_2  = __builtin_ia32_andpd(flag_2, y2);
#endif
#endif
  }




  typedef dvec real;
}

#define LDA(x)   (SIMD::load<true >(x))
#define LDU(x)   (SIMD::load<false>(x))
#define STU(x,y) (SIMD::store<false>(x, y))
#define STA(x)   (*(SIMD::real*)&x)

