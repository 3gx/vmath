#pragma once

#ifndef __SSE__
#error "Please enable SSE support"
#endif

#define __out

#include <cassert>

namespace SIMD
{
  enum {WIDTH = 2};
  typedef double _v2df  __attribute__((vector_size(16)));

  struct dvec 
  {
    _v2df val;
    dvec() {}
    dvec(_v2df _val) : val(_val) {};
    dvec(const double  f) : val((_v2df){f,f}) {}
    dvec(const double  a, const double b) : val((_v2df){a,b}) {}
    dvec(const double *p, const bool aligned = true) : 
      val(aligned ? *(_v2df*)p : __builtin_ia32_loadupd(p)) {}
    
    operator _v2df() const {return val;}

		const dvec& operator =(const dvec a) {val  = a.val; return *this;}
		const dvec& operator+=(const dvec a) {val += a.val; return *this;}
		const dvec& operator-=(const dvec a) {val -= a.val; return *this;}
		const dvec& operator*=(const dvec a) {val *= a.val; return *this;}
		const dvec& operator/=(const dvec a) {val /= a.val; return *this;}
		const dvec operator -()  const {return dvec(-val);}
		const dvec operator +(const dvec a) const {return dvec(val + a.val);}
		const dvec operator -(const dvec a) const {return dvec(val - a.val);}
		const dvec operator *(const dvec a) const {return dvec(val * a.val);}
		const dvec operator /(const dvec a) const {return dvec(val / a.val);}
		
    const dvec operator!=(const dvec a) const 
    {
			return dvec(__builtin_ia32_cmpneqpd(val, a.val));
		}
		const dvec operator==(const dvec a) const {
			return dvec(__builtin_ia32_cmpeqpd(val, a.val));
		}

		const dvec operator<(const dvec &rhs) const{
			return __builtin_ia32_cmpltpd(val, rhs.val);
		}
		const dvec operator<=(const dvec &rhs) const{
			return __builtin_ia32_cmplepd(val, rhs.val);
		}
		const dvec operator>(const dvec &rhs) const{
			return __builtin_ia32_cmpgtpd(val, rhs.val);
		}
		const dvec operator>=(const dvec &rhs) const{
			return __builtin_ia32_cmpgepd(val, rhs.val);
		}
		const dvec operator|(const dvec &rhs) const{
			return __builtin_ia32_orpd(val, rhs.val);
		}
		const dvec operator&(const dvec &rhs) const{
			return __builtin_ia32_andpd(val, rhs.val);
		}

    operator bool() const {return true;}
		
    const double operator[] (const int i) const {
			switch(i) {
				case 0:	return __builtin_ia32_vec_ext_v2df(val, 0);
				case 1:	return __builtin_ia32_vec_ext_v2df(val, 1);
				default: assert(0);
			}
    }
	

  };

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

  template<const bool ALIGN, typename T>
    inline static dvec store(T &ptr, const dvec val)
    {
      if (ALIGN)
      {
        *(_v2df*)&ptr = val;
      }
      else
      {
        __builtin_ia32_storeupd((double*)&ptr, val);
      }
      return val;
    }
  
  template<const bool ALIGN, typename T>
    inline static dvec load(const T &ptr)
    {
      return dvec((double*)&ptr, ALIGN);
    }

  template<const int N>
    inline static dvec broadcast(const dvec x)
    {
      const int mask =
        N == 0 ? (0 << 0) + (0 << 0) : (1 << 0) + (1 << 1);
      return __builtin_ia32_shufpd(x, x, mask);
    }



  typedef dvec real;
}

#define LDA(x)   (SIMD::load<true >(x))
#define LDU(x)   (SIMD::load<false>(x))
#define STU(x,y) (SIMD::store<false>(x, y))
#define STA(x)   (*(SIMD::real*)&x)

