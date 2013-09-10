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
  assert((nx % vreal::WIDTH*N) == 0);
  assert((ny % vreal::WIDTH*N) == 0);

  static VEC u[NX] __attribute__((aligned(64))) ;
  static VEC d[NTHREADS][N*NX] __attribute__((aligned(64))) ;

  Tridiag<VEC>::Reduce(nx, mat, u);

#pragma omp parallel
  {
    const int tid = omp_get_thread_num();

#pragma omp for nowait
    for (int kc = 0; kc < ny; kc += vreal::WIDTH*N)
    {
      for (int jc = 0; jc < nx; jc += vreal::WIDTH*N)
      {
        for (int k = 0; k < vreal::WIDTH*N; k++)
          for (int j = 0; j < N; j++)
            d[tid][jc*N+k*N+j] = VLDA(qin(jc+j*vreal::WIDTH, kc+k));
        SIMD::transpose<N>(&d[tid][jc*N]);
      }

      Tridiag<VEC>::template Solve<N>(nx, mat, u, d[tid]);


      for (int jc = 0; jc < nx; jc += vreal::WIDTH*N)
      {
        SIMD::transpose<N>(&d[tid][jc*N]);
        for (int k = 0; k < vreal::WIDTH*N; k++)
          for (int j = 0; j < N; j++)
            VREF(qout(jc+j*vreal::WIDTH, kc+k)) = d[tid][jc*N+k*N+j];
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
  const int NCACHE = MAX(4*vreal::WIDTH, NCACHEREAL);
  assert( NCACHE >= vreal::WIDTH);
  assert((NCACHE  % vreal::WIDTH) == 0);

  assert( nx > 0);
  assert( ny > 0);
  assert((nx % NCACHE)  == 0);

  const int NCACHESIMD = NCACHE/vreal::WIDTH;

  static vreal d[NTHREADS][NY*NCACHESIMD] __attribute__((aligned(64)));


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
          const vreal c = cnst[kc];
          for (int j = 0; j < NCACHESIMD; j++)
            d[tid][kc*NCACHESIMD+j] = VLDA(qin(jc+j*vreal::WIDTH, kc)) * c;
        }
      }
      else
      {
        for (int kc = 1; kc < ny-1; kc++)
        {
          const vreal c = cnst[kc];
          for (int j = 0; j < NCACHESIMD; j++)
            d[tid][kc*NCACHESIMD+j] = VLDA(qin(jc+j*vreal::WIDTH, kc)) * c;
        }
      }

      Tridiag<VEC>::template Solve<NCACHESIMD>(ny, mat, d[tid]);

      for (int kc = 0; kc < ny; kc++)
        for (int j = 0; j < NCACHESIMD; j++)
          VREF(qout(jc+j*vreal::WIDTH,kc)) += d[tid][kc*NCACHESIMD+j];
    }
  }

}
