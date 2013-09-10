#pragma once

/***********************************************************************
 *  this subroutine perform the calculation of dph , periodic direction
 *  along x3 and x1to use the real fourier transform
 ***********************************************************************/

void phcalcTP() /* check */
{
  static FFTW_COMPLEX xa[NTHREADS][NX];
  static real         xr[NTHREADS][NX];

  //   fft applied to the x2 direction to the
  //   complex coeff. from cos fft
  //   from physical to wave number space
  assert((n2m & 1) == 0);
  const int n2mh = n2m/2;
  const real cinv = 1.0/(real)n2m;

#pragma omp parallel
  {
    const int tid = omp_get_thread_num();

#pragma omp for
    for (int k = 0; k < n3m; k++)
    {
      for (int j = 0; j < n2m; j++)
        xr[tid][j] = dph(j,k);

      FFTW_EXEC_R2C(fwd_plan, xr[tid], xa[tid]);

      for (int j = 0; j < n2mh+1; j++)
        dph(j,k) = xa[tid][j][0] * cinv;

      for (int j = 0; j < n2mh+1; j++)
        dph(j+n2mh+1,k) = xa[tid][j][1] * cinv;
    }

    //m==================================================================
    //     inversion of the matrix in the x2 and x3 directions (real part)

    static Tridiag<real> mat[NTHREADS][NY];
    static real d[NTHREADS][NY];


#if 0
    for (int j = 0; j < n2m; j++)
    {
      const int jmh = j <= n2mh ? j : j - (n2mh+1);
      assert(jmh == jmhv(j)-1);

      for (int k = 0; k < n3m; k++)
      {
#if 0
        const real acphT_b = acphk(k) - ak2(jmh);
#else  /* egaburov: see my comment in phiniTP_intel.f:33 */
        const real acphT_b = acphk(k) - ak2(j);
#endif
        const real inv = 1.0/acphT_b;
        mat [k].a = amphk(k) * inv;
        mat [k].b = 1.0;
        mat [k].c = apphk(k) * inv;
        d   [k] = dph(j,k)   * inv;
      }

      Tridiag<real>::Reduce(n3m, mat);
      Tridiag<real>::Solve (n3m, mat, d);

      for (int k = 0; k < n3m; k++)
        dph(j,k) = d[k];
    }
#else
    {
      const int NCACHE = MAX(4*SIMD::WIDTH, NCACHEREAL);
      assert( NCACHE >= SIMD::WIDTH);
      assert((NCACHE  % SIMD::WIDTH) == 0);

      assert( n2m > 0);
      assert( n3m > 0);
      assert((n2m  % NCACHE)  == 0);

      const int NCACHESIMD = NCACHE/SIMD::WIDTH;

      static Tridiag<SIMD::real> mat[NTHREADS][NY*NCACHESIMD];
      static         SIMD::real    d[NTHREADS][NY*NCACHESIMD];

#pragma omp for
      for (int jc = 0; jc < n2m; jc += NCACHE)
      {
        for (int kc = 0; kc < n3m; kc++)
        {
          const SIMD::real a(amphk(kc));
          const SIMD::real b(acphk(kc));
          const SIMD::real c(apphk(kc));
          for (int j = 0; j < NCACHESIMD; j++)
          {
#if 0
            const SIMD::real inv = SIMD::real(1.0)/(b - LDA(ak2(jc+j*SIMD::WIDTH)));
            mat[tid][kc*NCACHESIMD+j].a = inv * a;
            mat[tid][kc*NCACHESIMD+j].b = 1.0    ;
            mat[tid][kc*NCACHESIMD+j].c = inv * c;
            d  [tid][kc*NCACHESIMD+j]   = inv * LDA(dph(jc+j*SIMD::WIDTH, kc));
#else
            mat[tid][kc*NCACHESIMD+j].a = a;
            mat[tid][kc*NCACHESIMD+j].b = b - LDA(ak2(jc+j*SIMD::WIDTH));
            mat[tid][kc*NCACHESIMD+j].c = c;
            d  [tid][kc*NCACHESIMD+j]   =     LDA(dph(jc+j*SIMD::WIDTH, kc));
#endif
          }
        }

        Tridiag<SIMD::real>::template ReduceN<NCACHESIMD>(n3m, mat[tid]);
        Tridiag<SIMD::real>::template  SolveN<NCACHESIMD>(n3m, mat[tid], d[tid]);

        for (int kc = 0; kc < n3m; kc++)
          for (int j = 0; j < NCACHESIMD; j++)
            STA(dph(jc+j*SIMD::WIDTH,kc)) = d[tid][kc*NCACHESIMD+j];
      }
    }

#pragma omp for
    for (int j = n2m; j < n2+1; j++)
    {
      const int jmh = jmhv(j) - 1;

      for (int k = 0; k < n3m; k++)
      {
#if 0
        const real inv = 1.0/(acphk(k) - ak2(jmh));
        mat[tid][k].a = amphk(k) * inv;
        mat[tid][k].b = 1.0;
        mat[tid][k].c = apphk(k) * inv;
        d  [tid][k]   = dph(j,k)   * inv;
#else
        mat[tid][k].a = amphk(k);
        mat[tid][k].b = acphk(k) - ak2(jmh);
        mat[tid][k].c = apphk(k);
        d  [tid][k]   = dph(j,k);
#endif
      }

      Tridiag<real>::Reduce(n3m, mat[tid]);
      Tridiag<real>::Solve (n3m, mat[tid], d[tid]);

      for (int k = 0; k < n3m; k++)
        dph(j,k) = d[tid][k];
    }
#endif

    //================================================================
    //   inverse fft applied to the phi x1 direction
    //   from wave number space to physical space

#pragma omp for
    for (int k = 0; k < n3m; k++)
    {
      for (int j = 0; j < n2mh+1; j++)
      {
        xa[tid][j][0] = dph(j,        k);
        xa[tid][j][1] = dph(j+n2mh+1, k);
      }

      FFTW_EXEC_C2R(bck_plan, xa[tid], xr[tid]);

      for (int j = 0; j < n2m; j++)
        dph(j,k) = xr[tid][j];
    }
  }

}
