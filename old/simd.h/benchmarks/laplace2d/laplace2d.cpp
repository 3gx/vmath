#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <cmath>

#include "defs.h"
#include "objs/laplace2d_ispc.h"

#include <sys/time.h>
static inline double rtc()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1.e-6 * tv.tv_usec;
}

template<typename T>
static inline T SQR(const T x) { return x*x; }

static int laplace2d(
    const int nx, const int ny,
    const real dataIn[], real dataOut[])
{
  const real (* in)[nx] = (real (*)[nx])&dataIn [0];
        real (*out)[nx] = (real (*)[nx])&dataOut[0];

  out[   0][   0] = in[   0][   0];
  out[   0][nx-1] = in[   0][nx-1];

#pragma omp parallel for schedule(dynamic)
  for (int j = 1; j < ny-1; j++)
    for (int i = 1; i < nx-1; i++)
      out[j][i] = (in[j][i-1] + in[j][i+1] + in[j-1][i] + in[j+1][i]) * real(0.25);

  out[ny-1][   0] = in[ny-1][   0];
  out[ny-1][nx-1] = in[ny-1][nx-1];
}

int main(int argc, char * argv [])
{
  const int nx = argc > 1 ? atoi(argv[1]) : 2048;
  const int ny = argc > 2 ? atoi(argv[2]) : 2048;

  fprintf(stderr, "laplace2d benchmark:  square domain %d x %d \n", nx, ny);

  const int size = nx*ny;

#if 0 
  std::vector<real> in_vec(size, 0.0), out_vec(size,0.0);
#else
  const int N = 2048*2049;
  assert(size < N);
  static real  in_vec[N] __attribute__((aligned(32)));
  static real out_vec[N] __attribute__((aligned(32)));
  static real   data0[N];
  static real   res_cpu[N];
#endif

  real (* in)[nx] = (real (*)[nx])& in_vec[0];
  real (*out)[nx] = (real (*)[nx])&out_vec[0];

  in[0][0] = in[0][nx-1] = in[ny-1][0] = in[ny-1][nx-1] = 1.0;

  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      in[j][i] = drand48();

  for (int i = 0; i < size; i++)
    data0[i] = in_vec[i];

  const int nRep = 10;
  {
    for (int i = 0; i < size; i++)
      in_vec[i]= data0[i];
    const double t0 = rtc();
    for (int r = 0; r < nRep; r++)
    {
      laplace2d(nx,ny,  &in_vec[0], &out_vec[0]);
      laplace2d(nx,ny, &out_vec[0],  &in_vec[0]);
    }
    const double t1 = rtc();

    double ckSum = 0.0;
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
        ckSum += in[j][i]*in[j][i];

    fprintf(stderr, " --> CPU :  done in %g sec :: BW= %g GB/s , ckSum = %g \n",
        t1-t0, nRep*2.0*size*sizeof(real)*2.0/(t1-t0)/1e9, ckSum/size);
  }

  {
    for (int i = 0; i < size; i++)
      in_vec[i]= data0[i];
    const double t0 = rtc();
    for (int r = 0; r < nRep; r++)
    {
      laplace2d(nx,ny,  &in_vec[0], &out_vec[0]);
      laplace2d(nx,ny, &out_vec[0],  &in_vec[0]);
    }
    const double t1 = rtc();

    double ckSum = 0.0;
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
        ckSum += in[j][i]*in[j][i];
    
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
      ispc::laplace2d(nx,ny,  &in_vec[0], &out_vec[0]);
      ispc::laplace2d(nx,ny, &out_vec[0],  &in_vec[0]);
    }
    const double t1 = rtc();

    double ckSum = 0.0;
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
        ckSum += in[j][i]*in[j][i];

    double err = 0.0;
    for (int i = 0; i < size; i++)
      err += SQR(in_vec[i] - res_cpu[i])/SQR(res_cpu[i]);


    fprintf(stderr, " --> ISPC:  done in %g sec :: BW= %g GB/s , ckSum = %g \n",
        t1-t0, nRep*2.0*size*sizeof(real)*2.0/(t1-t0)/1e9, ckSum/size);

    fprintf(stderr, " err= %g \n", std::sqrt(err));
  }


  return 0;
}
