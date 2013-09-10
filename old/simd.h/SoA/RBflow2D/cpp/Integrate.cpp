#include "Integrate.h"
#include "../dims.h"

#define CACHELINE (64*2)
#define NTHREADSMAX (16)


#ifdef SINGLE
typedef float real;
#else
typedef double real;
#endif

namespace ARR
{
  enum
  {
    RHS1 = 0,
    RHS2,
    RHS3,
    Q2, Q3, PR, DENS,
    HRO, RHS,
    RU1, RU2, RU3, RURO,
    DQ, QCAP, DPH,
    N
  };
}
static real tmpData[ARR::N][__M2__*__M3__] __attribute__((aligned(32)));

#define NCACHEREAL (CACHELINE/sizeof(real))


Problem::Integrate<real, __M2__, __M3__, NCACHEREAL, NTHREADSMAX> *prob = NULL;

extern "C"
{
  void lib_open_(
      real *dx2, real *dx3,
      real *rc, real *rm, real *g2rc, real *g2rm,
      real *zz, real *zm, real *g3rc, real *g3rm,

      real *udx2rm, real *udx2c, real *udx2m,
      real *usg2rc, real *udx3c, real *udx3m,

      int  *jmv, int *jmc, int *jpv, int *jpc, int *jmhv,
      int  *kmc, int *kpc, int *kmv, int *kpv, int *kup, int *kum,

      real *ap2j,  real *ac2j,  real *am2j,
      real *ap2je, real *ac2je, real *am2je,
      real *ap3j,  real *ac3j,  real *am3j,
      real *apscj, real *acscj, real *amscj,

      real *ap3ck,  real *ac3ck,  real *am3ck,
      real *ap3sk,  real *ac3sk,  real *am3sk,
      real *ap3ssk, real *ac3ssk, real *am3ssk,

      real *ak2,
      real *amphj, real *apphj, real *acphjj,
      real *amphk, real *acphk, real *apphk,

      int  *n2,  int  *n3,
      real *ren, real *pec,
      real *alm, real *gam, real *rom,
      real *rhs1, real *rhs2, real *rhs3)
      {
        assert(prob == NULL);

        int nthreads = 1;
#pragma omp parallel 
#pragma omp single
        nthreads = omp_get_num_threads();
        assert(nthreads <= NTHREADSMAX);

        prob = new Problem::Integrate<real, __M2__, __M3__, NCACHEREAL, NTHREADSMAX>(
            *dx2, *dx3,
            rc, rm, g2rc, g2rm,
            zz, zm, g3rc, g3rm,

            udx2rm, udx2c, udx2m,
            usg2rc, udx3c, udx3m,

            jmv, jmc, jpv, jpc, jmhv,
            kmc, kpc, kmv, kpv, kup, kum,

            ap2j,  ac2j,  am2j,
            ap2je, ac2je, am2je,
            ap3j,  ac3j,  am3j,
            apscj, acscj, amscj,

            ap3ck,  ac3ck,  am3ck,
            ap3sk,  ac3sk,  am3sk,
            ap3ssk, ac3ssk, am3ssk,

            ak2,
            amphj, apphj, acphjj,
            amphk, acphk, apphk,

            *n2, *n3,
            *ren, *pec,
            alm, gam, rom,

            tmpData[ARR::RHS1], tmpData[ARR::RHS2], tmpData[ARR::RHS3]
              );
#pragma omp parallel for
        for (int i = 0; i < __M2__*__M3__; i++)
        {
          tmpData[ARR::RU1 ][i] = 0.0;
          tmpData[ARR::RU2 ][i] = 0.0;
          tmpData[ARR::RU3 ][i] = 0.0;
          tmpData[ARR::RURO][i] = 0.0;
        }
      }

  void lib_close_()
  {
    assert(prob != NULL);
    delete prob;
    prob = NULL;
  }

  void lib_tschem_(real *dt, real *beta, int *nsst,
      real *denbs, real *denbn,
      real *q2, real *q3,
      real *pr, real *dens, real *hro, real *rhs,
      real *ru1, real *ru2,  real *ru3, real *ruro,
      real *dq, real *qcap,
      real *dph)
  {
    assert(prob != NULL);

#if 0  /* slower by 30%: ifort does something fishy with arrays on Westmere/SNB */
    prob->tschem(*dt, *beta, *nsst,
        denbs, denbn,
        q2, q3,
        pr, dens, hro, rhs,
        ru1, ru2, ru3, ruro,
        dq, qcap,
        dph);
#else

#if 0 /* copy overhead is too big */
#pragma omp parallel for
    for (int i = 0; i < __M2__*__M3__; i++)
    {
      tmpData[ARR::Q2  ][i] = q2  [i];
      tmpData[ARR::Q3  ][i] = q3  [i];
      tmpData[ARR::PR  ][i] = pr  [i];
      tmpData[ARR::DENS][i] = dens[i];
    }
    prob->tschem(*dt, *beta, *nsst,
        denbs, denbn,
        tmpData[ARR::Q2 ], tmpData[ARR::Q3  ], tmpData[ARR::PR ], tmpData[ARR::DENS],
        tmpData[ARR::HRO], tmpData[ARR::RHS ],
        tmpData[ARR::RU1], tmpData[ARR::RU2 ], tmpData[ARR::RU3], tmpData[ARR::RURO],
        tmpData[ARR::DQ ], tmpData[ARR::QCAP], tmpData[ARR::DPH]);
#pragma omp parallel for
    for (int i = 0; i < __M2__*__M3__; i++)
    {
      q2  [i] = tmpData[ARR::Q2  ][i];
      q3  [i] = tmpData[ARR::Q3  ][i];
      pr  [i] = tmpData[ARR::PR  ][i];
      dens[i] = tmpData[ARR::DENS][i];
    }
#else
    prob->tschem(*dt, *beta, *nsst,
        denbs, denbn,
        q2, q3, pr, dens,
        tmpData[ARR::HRO], tmpData[ARR::RHS ],
        tmpData[ARR::RU1], tmpData[ARR::RU2 ], tmpData[ARR::RU3], tmpData[ARR::RURO],
        tmpData[ARR::DQ ], tmpData[ARR::QCAP], tmpData[ARR::DPH]);
#endif

#endif

  }
};

