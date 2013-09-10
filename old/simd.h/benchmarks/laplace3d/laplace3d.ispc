#include "defs.h"


task void
laplace3d_task(
    uniform const int nytasks,
    uniform const int nztasks,
    uniform const int  nx,
    uniform const int  ny,
    uniform const int  nz,
    uniform const real dataIn[],
    uniform       real dataOut[])
{
#if 1
  const uniform int ytask = taskIndex % nytasks;
  const uniform int ztask = taskIndex / nytasks;
#else
  const uniform int ytask = taskIndex / nztasks;
  const uniform int ztask = taskIndex % nztasks;
#endif

  const uniform int ny_per_task = (ny-1) / nytasks + 1;
  const uniform int nyBeg = max(ytask * ny_per_task,    1);
  const uniform int nyEnd = min(nyBeg + ny_per_task, ny-1);
  
  const uniform int nz_per_task = (nz-1) / nztasks + 1;
  const uniform int nzBeg = max(ztask * nz_per_task,    1);
  const uniform int nzEnd = min(nzBeg + nz_per_task, nz-1);

#define  in(i,j,k)  dataIn[((k)*ny + (j))*nx + (i)]
#define out(i,j,k) dataOut[((k)*ny + (j))*nx + (i)]

  const uniform real a = 1.0;
  const uniform real b = 6.0;
  const uniform real sixth = a/b;

  for (uniform int k = nzBeg; k < nzEnd; k++)
    for (uniform int j = nyBeg; j < nyEnd; j++)
      foreach(i= 1 ... nx-1)
      {
#if 1
        out(i,j,k) = (
            in(i-1,j,k) + in(i+1,j,k) + 
            in(i,j-1,k) + in(i,j+1,k) + 
            in(i,j,k-1) + in(i,j,k+1) 
            ) * sixth;
#else
        out(i,j,k) = in(i,j,k);
#endif
      }

  if (taskIndex == 0)
  {
    out(   0,    0,    0) = in(   0,    0,    0);
    out(nx-1,    0,    0) = in(nx-1,    0,    0);
    out(   0, ny-1,    0) = in(   0, ny-1,    0);
    out(nx-1, ny-1,    0) = in(nx-1, ny-1,    0);
    out(   0,    0, nz-1) = in(   0,    0, nz-1);
    out(nx-1,    0, nz-1) = in(nx-1,    0, nz-1);
    out(   0, ny-1, nz-1) = in(   0, ny-1, nz-1);
    out(nx-1, ny-1, nz-1) = in(nx-1, ny-1, nz-1);
  }

#undef in
#undef out
}


export void laplace3d(
    uniform const int  nx,
    uniform const int  ny,
    uniform const int  nz,
    uniform       real dataIn[],
    uniform       real dataOut[])
{
  const uniform int nytasks = 8;
  const uniform int nztasks = 8;
  launch[nytasks * nztasks] laplace3d_task(
      nytasks, nztasks,
      nx, ny, nz,
      dataIn, dataOut);
}
