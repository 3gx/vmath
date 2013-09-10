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


#pragma omp barrier

    //m==================================================================
    //     inversion of the matrix in the x2 and x3 directions (real part)

    static Tridiag<real> mat[NTHREADS][NY];
    static real d[NTHREADS][NY];


    {
      const int NCACHE = MAX(4*vreal::WIDTH, NCACHEREAL);
      assert( NCACHE >= vreal::WIDTH);
      assert((NCACHE  % vreal::WIDTH) == 0);

      assert( n2m > 0);
      assert( n3m > 0);
      assert((n2m  % NCACHE)  == 0);

      const int NCACHESIMD = NCACHE/vreal::WIDTH;

      static Tridiag<vreal> mat[NTHREADS][NY*NCACHESIMD];
      static         vreal    d[NTHREADS][NY*NCACHESIMD];

#pragma omp for
      for (int jc = 0; jc < n2m; jc += NCACHE)
      {
        for (int kc = 0; kc < n3m; kc++)
        {
          const vreal a(amphk(kc));
          const vreal b(acphk(kc));
          const vreal c(apphk(kc));
          for (int j = 0; j < NCACHESIMD; j++)
          {
#if 0
            const vreal inv = vreal(1.0)/(b - VLDA(ak2(jc+j*vreal::WIDTH)));
            mat[tid][kc*NCACHESIMD+j].a = inv * a;
            mat[tid][kc*NCACHESIMD+j].b = 1.0    ;
            mat[tid][kc*NCACHESIMD+j].c = inv * c;
            d  [tid][kc*NCACHESIMD+j]   = inv * VLDA(dph(jc+j*vreal::WIDTH, kc));
#else
            mat[tid][kc*NCACHESIMD+j].a = a;
            mat[tid][kc*NCACHESIMD+j].b = b - VLDA(ak2(jc+j*vreal::WIDTH));
            mat[tid][kc*NCACHESIMD+j].c = c;
            d  [tid][kc*NCACHESIMD+j]   =     VLDA(dph(jc+j*vreal::WIDTH, kc));
#endif
          }
        }

        Tridiag<vreal>::template ReduceN<NCACHESIMD>(n3m, mat[tid]);
        Tridiag<vreal>::template  SolveN<NCACHESIMD>(n3m, mat[tid], d[tid]);

        for (int kc = 0; kc < n3m; kc++)
          for (int j = 0; j < NCACHESIMD; j++)
            VREF(dph(jc+j*vreal::WIDTH,kc)) = d[tid][kc*NCACHESIMD+j];
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

#pragma omp barrier

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
