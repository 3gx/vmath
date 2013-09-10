#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>
#include <omp.h>

#include "defs.h"
#include "objs/laplace3d_ispc.h"

#include <sys/time.h>
static inline double rtc()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1.e-6 * tv.tv_usec;
}

template<typename T>
static inline T SQR(const T x) { return x*x; }

int threadId;
#pragma omp threadprivate(threadId)


static int laplace3d(
    const int nx, const int ny, const int nz,
    const real dataIn[], real dataOut[])
{
  const real (* in)[ny][nx] = (real (*)[ny][nx])&dataIn [0];
        real (*out)[ny][nx] = (real (*)[ny][nx])&dataOut[0];

#pragma omp parallel for schedule(static)
  for (int k = 1; k < nz-1; k++)
    for (int j = 1; j < ny-1; j++)
      for (int i = 1; i < nx-1; i++)
      {
#if 1
        out[k][j][i] = (
            in[k][j][i-1] + in[k][j][i+1] + 
            in[k][j-1][i] + in[k][j+1][i] + 
            in[k-1][j][i] + in[k+1][j][i] 
            ) * real(1.0/6.0);
#else
        out[k][j][i] = in[k][j][i];
#endif
      }

  out[   0][   0][   0] = in[   0][   0][   0];
  out[   0][nx-1][   0] = in[   0][nx-1][   0];
  out[ny-1][   0][   0] = in[ny-1][   0][   0];
  out[ny-1][nx-1][   0] = in[ny-1][nx-1][   0];
  out[   0][   0][nz-1] = in[   0][   0][nz-1];
  out[   0][nx-1][nz-1] = in[   0][nx-1][nz-1];
  out[ny-1][   0][nz-1] = in[ny-1][   0][nz-1];
  out[ny-1][nx-1][nz-1] = in[ny-1][nx-1][nz-1];
}

int main(int argc, char * argv [])
{
  const int nx = argc > 1 ? atoi(argv[1]) : 256;
  const int ny = argc > 2 ? atoi(argv[2]) : 256;
  const int nz = argc > 3 ? atoi(argv[3]) : 256;

  fprintf(stderr, "laplace3d benchmark:  square domain %d x %d x %d \n", nx, ny, nz);

  const int size = nx*ny*nz;

#if 0
  std::vector<real> in_vec(size, 0.0), out_vec(size,0.0), data0(size);
  std::vector<real> res_cpu(size);
#else
  real * in_vec = (real*)malloc(sizeof(real)*size);
  real *out_vec = (real*)malloc(sizeof(real)*size);
  real *data0   = (real*)malloc(sizeof(real)*size);
  real *res_cpu = (real*)malloc(sizeof(real)*size);
#endif

  int nthreads = 0;
#pragma omp parallel
#pragma omp critical
  nthreads++;

#pragma omp parallel
  threadId = omp_get_thread_num();


  fprintf(stderr, " number of OpenMP threads= %d\n", nthreads);


  real (* in)[ny][nx] = (real (*)[ny][nx])& in_vec[0];
  real (*out)[ny][nx] = (real (*)[ny][nx])&out_vec[0];


#if 1
#pragma omp parallel for schedule(static)
  for (int k = 1; k < nz-1; k++)
    for (int j = 1; j < ny-1; j++)
      for (int i = 1; i < nx-1; i++)
        in[k][j][i] = out[k][j][i] = 0.0;
#else
  //#pragma omp parallel for schedule(dynamic)
  for (int k = 1; k < nz-1; k++)
    for (int j = 1; j < ny-1; j++)
      for (int i = 1; i < nx-1; i++)
        in[k][j][i] = out[k][j][i] = 0.0;
#endif

#if 1  
#pragma omp parallel for schedule(static)
  for (int k = 0; k < nz; k++)
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
      {
        in [k][j][i] = drand48();
        out[k][j][i] = 0.0;
      }
#endif

  for (int i = 0; i < size; i++)
    data0[i] = in_vec[i];

  const int nRep = 10;
  {
#pragma omp parallel for
    for (int i = 0; i < size; i++)
      in_vec[i]= data0[i];

    const double t0 = rtc();
    for (int r = 0; r < nRep; r++)
    {
      laplace3d(nx,ny,nz,  &in_vec[0], &out_vec[0]);
      laplace3d(nx,ny,nz, &out_vec[0],  &in_vec[0]);
    }
    const double t1 = rtc();

    double ckSum = 0.0;
    for (int k = 0; k < ny; k++)
      for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
          ckSum += in[k][j][i]*in[k][j][i];

    fprintf(stderr, " --> CPU :  done in %g sec :: BW= %g GB/s , ckSum = %g \n",
        t1-t0, nRep*2.0*size*sizeof(real)*2.0/(t1-t0)/1e9, ckSum/size);
  }

  {
    for (int i = 0; i < size; i++)
      in_vec[i]= data0[i];
    const double t0 = rtc();
    for (int r = 0; r < nRep; r++)
    {
      laplace3d(nx,ny,nz,  &in_vec[0], &out_vec[0]);
      laplace3d(nx,ny,nz, &out_vec[0],  &in_vec[0]);
    }
    const double t1 = rtc();

    double ckSum = 0.0;
    for (int k = 0; k < ny; k++)
      for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
          ckSum += in[k][j][i]*in[k][j][i];

    for (int i = 0; i < size; i++)
      res_cpu[i] = in_vec[i];

    fprintf(stderr, " --> CPU :  done in %g sec :: BW= %g GB/s , ckSum = %g \n",
        t1-t0, nRep*2.0*size*sizeof(real)*2.0/(t1-t0)/1e9, ckSum/size);
  }

  {
    for (int i = 0; i < size; i++)
      in_vec[i]= data0[i];

    const double t0 = rtc();
    for (int r = 0; r < nRep; r++)
    {
      ispc::laplace3d(nx,ny,nz,  &in_vec[0], &out_vec[0]);
      ispc::laplace3d(nx,ny,nz, &out_vec[0],  &in_vec[0]);
    }
    const double t1 = rtc();

    double ckSum = 0.0;
    for (int k = 0; k < ny; k++)
      for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++)
          ckSum += in[k][j][i]*in[k][j][i];

    double err = 0.0;
    for (int i = 0; i < size; i++)
      err += SQR(in_vec[i] - res_cpu[i])/SQR(res_cpu[i]);



    fprintf(stderr, " --> ISPC:  done in %g sec :: BW= %g GB/s , ckSum = %g \n",
        t1-t0, nRep*2.0*size*sizeof(real)*2.0/(t1-t0)/1e9, ckSum/size);

    fprintf(stderr, " err= %g \n", std::sqrt(err));
  }


  return 0;
}
