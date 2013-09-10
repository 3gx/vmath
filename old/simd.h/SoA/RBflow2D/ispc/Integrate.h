#pragma once

#include "objs/init_ispc.h"
#include "objs/hdnl_ispc.h"
#include "objs/invtr_ispc.h"

#ifdef __OMP__
#include <omp.h>
#else
#define omp_get_thread_num() (0)
#define omp_get_num_threads() (1)
#endif

#include "FArray.h"
#include <fftw3.h>
#include <iostream>

#if 0
#define __mySSE__
#endif

#define MAX(x,y) (x>y ? x : y)

#if !defined(__mySSE__) && !defined(__AVX__)
#define __mySSE__
#endif


#ifdef __mySSE__

#ifdef SINGLE
#include "sse_fp32.h"
#else
#include "sse_fp64.h"
#endif

#else  /* AVX */

#ifdef SINGLE
#include "avx_fp32.h"
#else
#include "avx_fp64.h"
#endif

#endif /* mySSE */

extern "C"
{
  double rtc(void);
}

#ifdef __PAPI__
#include <papi.h>
#define _PAPI_
#define START_PAPI {\
  start_usecs = PAPI_get_virt_usec(); \
  assert(PAPI_start_counters(Events, NUM_EVENTS) == PAPI_OK); }

#define STOP_PAPI {\
  assert(PAPI_stop_counters(values, NUM_EVENTS) == PAPI_OK); \
  stop_usecs = PAPI_get_virt_usec();}
#define TIME_PAPI(routine) {\
  START_PAPI \
  routine; \
  STOP_PAPI }
#else
#define START_PAPI
#define STOP_PAPI
#define TIME_PAPI
#endif

#ifdef SINGLE
#define FFTW_PLAN     fftwf_plan
#define FFTW_COMPLEX  fftwf_complex
#define FFTW_PLAN_R2C fftwf_plan_dft_r2c_1d
#define FFTW_PLAN_C2R fftwf_plan_dft_c2r_1d
#define FFTW_EXEC_R2C fftwf_execute_dft_r2c
#define FFTW_EXEC_C2R fftwf_execute_dft_c2r
#else
#define FFTW_PLAN     fftw_plan
#define FFTW_COMPLEX  fftw_complex
#define FFTW_PLAN_R2C fftw_plan_dft_r2c_1d
#define FFTW_PLAN_C2R fftw_plan_dft_c2r_1d
#define FFTW_EXEC_R2C fftw_execute_dft_r2c
#define FFTW_EXEC_C2R fftw_execute_dft_c2r
#endif

#if 0
#define _GFORTRAN_
#endif

#ifdef _GFORTRAN_
extern double __param_MOD_mytimer[];
#define mytimer __param_MOD_mytimer
#else /* _IFORT_ */
extern double param_mp_mytimer_[];
#define mytimer param_mp_mytimer_
#endif

template<class T>
inline T __min(const T a, const T b) {return std::min(a,b);}
template<class T>
inline T __max(const T a, const T b) {return std::max(a,b);}

#define CCLOCK (1.0e9)
inline static unsigned long long TIMER()
{
  unsigned int low,high;
  __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high));
  unsigned long long count = (unsigned long long)low + (((unsigned long long) high)<<32);
  return count;
}

namespace Problem
{
  template < int start, int end, typename T >                                                          
    struct Loop                                                                                          
    {                                                                                                    
      static void eval()     {                                                                         
        const T v(start);                    // do work                                              
        Loop<start+1,end, T>::eval() ;   // increment argument and                 // recurse        
      }                                                                                                
    };                                                                                                   
  template <int end, typename T >                                                                      
    struct Loop <end, end, T>                                                                            
    {                                                                                                    
      static void eval() {}
    };        

  template<typename real, const int NX, const int NY, const int NCACHEREAL, const int NTHREADS>
    struct Integrate
    {
      FArray2D<real, NX,    NY> rhs1, rhs2, rhs3;

      /* local array */

      FArray2D<real, NX,    NY> q2, q3;
      FArray2D<real, NX,    NY> pr,  dens, hro, rhs;
      FArray2D<real, NX,    NY> ru1, ru2,  ru3, ruro;
      FArray2D<real, NX,    NY> dq, qcap;
      FArray2D<real, NX+16, NY> dph;

      /* grid parameters */

      const real dx2,  dx3;
      const real dx2q, dx3q;
      const FArray1D<real, NX> rc, rm, g2rc, g2rm;
      const FArray1D<real, NY> zz, zm, g3rc, g3rm;

      /******* QUANTITIES FOR DERIVATIVES ********/

      const FArray1D<real, NX> udx2rm, udx2c, udx2m;
      const FArray1D<real, NX> usg2rc;
      const FArray1D<real, NY> udx3c, udx3m;

      /******* Grid indices**************************************/

      const FArray1D<int, NX>   jmv,jmc,jpv,jpc;
      const FArray1D<int, NX+1> jmhv;
      const FArray1D<int, NY>   kmc,kpc,kmv,kpv,kup,kum;

      /******* Metric coefficients *******************************/

      const FArray1D<real, NX> ap2j,ac2j,am2j;
      const FArray1D<real, NX> ap2je,ac2je,am2je;
      const FArray1D<real, NX> ap3j,ac3j,am3j;
      const FArray1D<real, NX> apscj,acscj,amscj;

      const FArray1D<real, NY> ap3ck,ac3ck,am3ck;
      const FArray1D<real, NY> ap3sk,ac3sk,am3sk;
      const FArray1D<real, NY> ap3ssk,ac3ssk,am3ssk;

      /******* Variables for FFT and Poisson solver****************/

      const FArray1D<real, NX> ak2;
      const FArray1D<real, NX> amphj,apphj,acphjj;
      const FArray1D<real, NY> amphk,apphk,acphk;

      FFTW_COMPLEX xa_vec[NX];
      real         xr_vec[NX];

      FArray1D <FFTW_COMPLEX, NX> xa_arr;
      FArray1D <real,         NX> xr_arr;

      FFTW_PLAN fwd_plan, bck_plan;

      /******* Other variables ***********************************/

      const int n2, n3;
      const int n2m, n3m;
      const real ren, pec;

      const FArray1D<real,3> alm, gam, rom;
      real al,ga,ro;
      FArray1D<real, NX> denbs, denbn;

      real dt, beta;

      template<typename T>
        inline static T SQR(const T x) {return x*x;}

#include "tridiag.h"  /* done */
#include "solT.h"

#include "hdnl.h"     /* done */
#include "hdnl_all.h" /* done */
#include "invtr2.h"   /* done */
#include "solq2.h"    /* done */
#include "invtr3.h"   /* done */
#include "solq3.h"    /* done */
#include "invtrro.h"  /* done */
#include "solro.h"    /* done */
#include "divg.h"     /* done */
#include "updvp.h"    /* done */
#include "prcalc.h"   /* done */
#include "phcalcTP.h" /* done */
#include "invtr_all.h"   /* done */
#include "update_qp.h"


      Integrate(
          const real _dx2, const real _dx3,
          const real __restrict__ *_rc, const real __restrict__ *_rm, const real __restrict__ *_g2rc, const real __restrict__ *_g2rm,
          const real __restrict__ *_zz, const real __restrict__ *_zm, const real __restrict__ *_g3rc, const real __restrict__ *_g3rm,

          const real __restrict__ *_udx2rm, const real __restrict__ *_udx2c, const real __restrict__ *_udx2m,
          const real __restrict__ *_usg2rc, const real __restrict__ *_udx3c, const real __restrict__ *_udx3m,

          const int   *_jmv, const int  *_jmc, const int  *_jpv, const int  *_jpc, const int  *_jmhv,
          const int   *_kmc, const int  *_kpc, const int  *_kmv, const int  *_kpv, const int  *_kup, const int  *_kum,

          const real __restrict__ *_ap2j,  const real __restrict__ *_ac2j,  const real __restrict__ *_am2j,
          const real __restrict__ *_ap2je, const real __restrict__ *_ac2je, const real __restrict__ *_am2je,
          const real __restrict__ *_ap3j,  const real __restrict__ *_ac3j,  const real __restrict__ *_am3j,
          const real __restrict__ *_apscj, const real __restrict__ *_acscj, const real __restrict__ *_amscj,

          const real __restrict__ *_ap3ck,  const real __restrict__ *_ac3ck,  const real __restrict__ *_am3ck,
          const real __restrict__ *_ap3sk,  const real __restrict__ *_ac3sk,  const real __restrict__ *_am3sk,
          const real __restrict__ *_ap3ssk, const real __restrict__ *_ac3ssk, const real __restrict__ *_am3ssk,

          const real __restrict__ *_ak2,
          const real __restrict__ *_amphj, const real __restrict__ *_apphj, const real __restrict__ *_acphjj,
          const real __restrict__ *_amphk, const real __restrict__ *_acphk, const real __restrict__ *_apphk,

          const int  _n2,  const int  _n3,
          const real _ren, const real _pec, 
          const real __restrict__ *_alm,  const real __restrict__ *_gam, const real __restrict__ *_rom,
          real *rhs1a, real *rhs2a, real *rhs3a)
            :
              rhs1(rhs1a), rhs2(rhs2a), rhs3(rhs3a),
              dx2(_dx2), dx3(_dx3), dx2q(dx2*dx2), dx3q(dx3*dx3),
              rc(_rc), rm(_rm), g2rc(_g2rc), g2rm(_g2rm),
              zz(_zz), zm(_zm), g3rc(_g3rc), g3rm(_g3rm),

              udx2rm(_udx2rm), udx2c(_udx2c), udx2m(_udx2m),
              usg2rc(_usg2rc), udx3c(_udx3c), udx3m(_udx3m),

              jmv(_jmv), jmc(_jmc), jpv(_jpv), jpc(_jpc), jmhv(_jmhv),
              kmc(_jmc), kpc(_jpc), kmv(_kmv), kpv(_kpv), kup(_kup), kum(_kum),

              ap2j (_ap2j ), ac2j (_ac2j ), am2j (_am2j ),
              ap2je(_ap2je), ac2je(_ac2je), am2je(_am2je),
              ap3j (_ap3j ), ac3j (_ac3j ), am3j (_am3j ),
              apscj(_apscj), acscj(_acscj), amscj(_amscj),

              ap3ck (_ap3ck ), ac3ck (_ac3ck ), am3ck (_am3ck ),
              ap3sk (_ap3sk ), ac3sk (_ac3sk ), am3sk (_am3sk ),
              ap3ssk(_ap3ssk), ac3ssk(_ac3ssk), am3ssk(_am3ssk),

              ak2(_ak2),
              amphj(_amphj), apphj(_apphj), acphjj(_acphjj),
              amphk(_amphk), apphk(_apphk), acphk (_acphk ),

              n2 (_n2 ), n3 (_n3 ), n2m(n2-1), n3m(n3-1),
              ren(_ren), pec(_pec),
              alm(_alm), gam(_gam), rom(_rom)
              {
#if 0
                assert(n2 == NX);
                assert(n3 == NY);
#endif
                xa_arr.set_data(xa_vec);
                xr_arr.set_data(xr_vec);

                assert(((n2-1)&1) == 0);
                fwd_plan = FFTW_PLAN_R2C(n2-1, xr_arr, xa_arr, FFTW_EXHAUSTIVE);
                bck_plan = FFTW_PLAN_C2R(n2-1, xa_arr, xr_arr, FFTW_EXHAUSTIVE);

                ispc::init(
                    (real*)_rc,  (real*)_rm,  (real*)_g2rc,  (real*)_g2rm,
                    (real*)_zz,  (real*)_zm,  (real*)_g3rc,  (real*)_g3rm,

                    (real*)_udx2rm,  (real*)_udx2c,  (real*)_udx2m,
                    (real*)_usg2rc,  (real*)_udx3c,  (real*)_udx3m,

                    (int*)_jmv, (int*)_jmc, (int*)_jpv, (int*)_jpc, (int*)_jmhv,
                    (int*)_kmc, (int*)_kpc, (int*)_kmv, (int*)_kpv, (int*)_kup, (int*)_kum,

                    (real*)_ap2j,   (real*)_ac2j,   (real*)_am2j,
                    (real*)_ap2je,  (real*)_ac2je,  (real*)_am2je,
                    (real*)_ap3j,   (real*)_ac3j,   (real*)_am3j,
                    (real*)_apscj,  (real*)_acscj,  (real*)_amscj,

                    (real*)_ap3ck,   (real*)_ac3ck,   (real*)_am3ck,
                    (real*)_ap3sk,   (real*)_ac3sk,   (real*)_am3sk,
                    (real*)_ap3ssk,  (real*)_ac3ssk,  (real*)_am3ssk,

                    (real*)_ak2,
                    (real*)_amphj,  (real*)_apphj,  (real*)_acphjj,
                    (real*)_amphk,  (real*)_acphk,  (real*)_apphk,

                    n2m,  n3m,
                    ren,  pec);
              }

      enum {NUM_EVENTS=3};
      long long values[NUM_EVENTS];
      int Events[NUM_EVENTS];

      long long values_long[NUM_EVENTS];
      long long time_long;
      unsigned long long start_usecs, stop_usecs;


#if 0
#define DEBUG
#endif

      template<const int NN, const int MM>
        void ckSum(FArray2D<real, NN, MM> &q, const char line[], const int inc = 0)
        {
#ifdef DEBUG
          double cksum = 0.0;
          for (int j = inc; j < n2m-inc; j++)
            for (int k = inc; k < n3m-inc; k++)
              cksum += q(j,k)*q(j,k);
          fprintf(stderr, "%s= %16.15g\n", line, cksum);
#endif
        }

      void tschem(const real dt, const real beta, const int nsst,
          real *denbs, real *denbn,
          real *q2,  real *q3,
          real *pr,  real *dens, real *hro, real *rhs,
          real *ru1, real *ru2,  real *ru3, real *ruro,
          real *dq,  real *qcap,
          real *dph)
      {


#ifdef _PAPI_
        static bool PAPIinit = false;
        static int cnt = 0;
        if (!PAPIinit)
        {
          Events[0] = PAPI_FP_OPS;
          Events[1] = PAPI_TOT_INS;
          Events[2] = PAPI_TOT_CYC;

          time_long = 0;
          for (int i = 0; i < NUM_EVENTS; i++)
            values_long[i] = 0;

          assert(PAPI_library_init(PAPI_VER_CURRENT) == PAPI_VER_CURRENT);
          start_usecs = stop_usecs = 0;

          PAPIinit = true;
        }
#endif

        this->dt   = dt;
        this->beta = beta;
        this->denbs.set_data(denbs);
        this->denbn.set_data(denbn);
        this->q2   .set_data(q2);
        this->q3   .set_data(q3);
        this->pr   .set_data(pr);
        this->dens .set_data(dens);
        this->hro  .set_data(hro);
        this->rhs  .set_data(rhs);
        this->ru1  .set_data(ru1);
        this->ru2  .set_data(ru2);
        this->ru3  .set_data(ru3);
        this->ruro .set_data(ruro);
        this->dq   .set_data(dq);
        this->qcap .set_data(qcap);
        this->dph  .set_data(dph);

        ispc::set_data(
            dt, beta,
            denbs, denbn,
            q2, q3,
            pr, dens, hro,
            rhs1, rhs2, rhs3,
            ru1, ru2, ru3, ruro,
            qcap, dph);



        static int val = 0;

        val = val + 1;


        double t0, t1;
        const int inc = 1;
        for (int ns = 0; ns < nsst; ns++)
        {
          al = alm(ns);
          ga = gam(ns);
          ro = rom(ns);

          ispc::set_algaro(al, ga, ro);


#ifdef DEBUG
          fprintf(stderr ,"  %d 0--------------------------0\n", val);
          ckSum(this->q2, "q2", inc); ckSum(this->q3, "q3", inc);
#endif

#if 1
          START_PAPI;
#endif

          t0 = rtc();

#if 0
          hdnl2();  /* vec, omp */
          t1 = rtc();  mytimer[0] += t1 - t0;  t0 = t1; ckSum(this->dph, "hdnl2::dph", inc);
          hdnl3();  /* vec, omp */
          t1 = rtc();  mytimer[1] += t1 - t0;  t0 = t1; ckSum(this->qcap, "hdnl3::qcap", inc);
          hdnlro();  /* vec, omp */
          t1 = rtc(); mytimer[2] += t1 - t0 ;  t0 = t1; ckSum(this->hro, "hdnlro::hro", inc);
#else  /*faster even in serial mode */
#if 0
          hdnl_all();  /* vec, omp */
#else
          ispc::hdnl_all();
#endif
          t1 = rtc();  mytimer[0] += t1 - t0;  t0 = t1; 
          ckSum(this->dph,  "hdnl_all::dph",  inc);
          ckSum(this->qcap, "hdnl_all::qcap", inc);
          ckSum(this->hro,  "hdnl_all::hro",  inc);
#endif

#if 0
          invtr2();   /* vec, omp */
          t1 = rtc(); mytimer[3] += t1 - t0;  t0 = t1; ckSum(this->q2, "invtr2::q2", inc);
          invtr3();   /* vec, omp */
          t1 = rtc();  mytimer[4] += t1 - t0; t0 = t1; ckSum(this->q3, "invtr3::q3", inc);
          invtrro();   /* vec, omp */
          t1 = rtc();  mytimer[8] += t1 - t0; t0 = t1; ckSum(this->dens, "invtrro::dens", inc);
#else
#if 0
          invtr_all(); /* vec, omp */
#else
          ispc::invtr_all();  /* some problem in this one ... */
#endif
          t1 = rtc(); mytimer[3] += t1 - t0;  t0 = t1; 

          double tA = rtc();
          solq2j();
          mytimer[49] += rtc() - tA;

          tA = rtc();
          solq2k();
          mytimer[50] += rtc() - tA;

          solq3j();
          solq3k();
          solroj();
          solrok();

          t1 = rtc(); mytimer[4] += t1 - t0;  t0 = t1; 
          ckSum(this->q2, "invtr2::q2  ", inc);
          ckSum(this->q3, "invtr2::q3  ", inc);
          ckSum(this->dens, "invtr2::dens", inc);
#endif

          divg();  /* vec, omp */
          t1 = rtc();  mytimer[5] += t1 - t0; t0 = t1; ckSum(this->dph, "divg::dph", inc);

          phcalcTP();   /* vec, omp */
          t1 = rtc(); mytimer[6] += t1 - t0;  t0 = t1; ckSum(this->dph, "phcalc::dph", inc);

#if 0
          updvp();  /* vec, omp */
          t1 = rtc(); mytimer[7] += t1 - t0;  t0 = t1; ckSum(this->q2, "updvp::q2", inc); ckSum(this->q3, "updvp::q3", inc);

          prcalc();  /* vec, omp */
          t1 = rtc();  mytimer[9] += t1 - t0; t0 = t1; ckSum(this->pr, "prcalc::pr", inc);
#else
          update_qp(); /* vec, omp */
          t1 = rtc(); mytimer[7] += t1 - t0;  t0 = t1; 
          ckSum(this->q2, "updvp::q2", inc); 
          ckSum(this->q3, "updvp::q3", inc);
          ckSum(this->pr, "updvp::pr", inc);
#endif

#if 1
          STOP_PAPI;
#endif

        }

#ifdef _PAPI_
        for (int i = 0; i < NUM_EVENTS; i++)
          values_long[i] += values[i];
        time_long += stop_usecs - start_usecs;
        if (++cnt == 10)
        {
          const double c = 1.0e6/1.0e9/(double)time_long;
          cnt = 0;
          fprintf(stderr, " >>> GFLOP/s= %g  GINS/s= %g  GCYC/s= %g : FPC= %g IPC = %g\n", 
              c*(double)values_long[0], c*(double)values_long[1], c*(double)values_long[2],
              (double)values_long[0]/(double)values_long[2],
              (double)values_long[1]/(double)values_long[2]);
        }
#endif
      }

    };
}
