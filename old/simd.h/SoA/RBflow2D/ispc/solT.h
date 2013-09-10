#pragma once

template<typename VEC>
void static solveX(
    const int nx, 
    const int ny,
    const FArray2D<real, NX, NY> &qin, 
    __out FArray2D<real, NX, NY> &qout,
    __out Tridiag<VEC> mat[]) 
{
  const int N = 2;
  assert( nx > 0);
  assert( ny > 0);
  assert((nx % SIMD::WIDTH*N) == 0);
  assert((ny % SIMD::WIDTH*N) == 0);

  static VEC u[NX] __attribute__((aligned(64))) ;
  static VEC d[NTHREADS][N*NX] __attribute__((aligned(64))) ;

  Tridiag<VEC>::Reduce(nx, mat, u);

#pragma omp parallel
  {
    const int tid = omp_get_thread_num();

#pragma omp for nowait
    for (int kc = 0; kc < ny; kc += SIMD::WIDTH*N)
    {
      for (int jc = 0; jc < nx; jc += SIMD::WIDTH*N)
      {
        for (int k = 0; k < SIMD::WIDTH*N; k++)
          for (int j = 0; j < N; j++)
            d[tid][jc*N+k*N+j] = SIMD::real(&qin(jc+j*SIMD::WIDTH, kc+k), true);
        SIMD::transpose<N>(&d[tid][jc*N]);
      }

      Tridiag<VEC>::template Solve<N>(nx, mat, u, d[tid]);


      for (int jc = 0; jc < nx; jc += SIMD::WIDTH*N)
      {
        SIMD::transpose<N>(&d[tid][jc*N]);
        for (int k = 0; k < SIMD::WIDTH*N; k++)
          for (int j = 0; j < N; j++)
            *(SIMD::real*)&qout(jc+j*SIMD::WIDTH, kc+k) = d[tid][jc*N+k*N+j];
      }
    }
  }
}


template<typename VEC>
static void solveY(
    const int nx, 
    const int ny,
    const VEC cnst[],
    const FArray2D<real, NX, NY> &qin, 
    __out FArray2D<real, NX, NY> &qout,
    __out Tridiag<VEC> mat[],
    const bool flag = false) 
{
  const int NCACHE = MAX(4*SIMD::WIDTH, NCACHEREAL);
  assert( NCACHE >= SIMD::WIDTH);
  assert((NCACHE  % SIMD::WIDTH) == 0);

  assert( nx > 0);
  assert( ny > 0);
  assert((nx % NCACHE)  == 0);

  const int NCACHESIMD = NCACHE/SIMD::WIDTH;

  static SIMD::real d[NTHREADS][NY*NCACHESIMD] __attribute__((aligned(64)));


  Tridiag<VEC>::Reduce(nx, mat);
#pragma omp parallel
  {
    const int tid = omp_get_thread_num();

    if (flag)
      for (int j = 0; j < NCACHESIMD; j++)
        d[tid][0+j] = d[tid][(ny-1)*NCACHESIMD+j] = 0.0;

#pragma omp for 
    for (int jc = 0; jc < nx; jc += NCACHE)
    {
      if (!flag)
      {
        for (int kc = 0; kc < ny; kc++)
        {
          const SIMD::real c = cnst[kc];
          for (int j = 0; j < NCACHESIMD; j++)
            d[tid][kc*NCACHESIMD+j] = SIMD::real(&qin(jc+j*SIMD::WIDTH, kc), true) * c;
        }
      }
      else
      {
        for (int kc = 1; kc < ny-1; kc++)
        {
          const SIMD::real c = cnst[kc];
          for (int j = 0; j < NCACHESIMD; j++)
            d[tid][kc*NCACHESIMD+j] = SIMD::real(&qin(jc+j*SIMD::WIDTH, kc), true) * c;
        }
      }

      Tridiag<VEC>::template Solve<NCACHESIMD>(ny, mat, d[tid]);

      for (int kc = 0; kc < ny; kc++)
        for (int j = 0; j < NCACHESIMD; j++)
          *(SIMD::real*)&qout(jc+j*SIMD::WIDTH,kc) += d[tid][kc*NCACHESIMD+j];
    }
  }

}
