#include "defs.h"


task void
laplace2d_task(
    uniform const int ntasks,
    uniform const int  nx,
    uniform const int  ny,
    uniform const real dataIn[],
    uniform       real dataOut[])
{
  const uniform int ny_per_task = (ny-1)/ntasks + 1;
  const uniform int nyBeg = max(taskIndex * ny_per_task,   1);
  const uniform int nyEnd = min(nyBeg + ny_per_task,    ny-1);

#define  in(i,j)  dataIn[(j)*nx + (i)]
#define out(i,j) dataOut[(j)*nx + (i)]

  for (uniform int j = nyBeg; j < nyEnd; j++)
    foreach(i= 1 ... nx-1)
      out(i,j) = (in(i-1,j) + in(i+1,j) + in(i,j-1) + in(i,j+1)) * (real)0.25;

  if (taskIndex == 0)
  {
    out(   0,    0) = in(   0,    0);
    out(nx-1,    0) = in(nx-1,    0);
    out(   0, ny-1) = in(   0, ny-1);
    out(nx-1, ny-1) = in(nx-1, ny-1);
  }

#undef in
#undef out
}

export void laplace2d(
    uniform const int  nx,
    uniform const int  ny,
    uniform const real dataIn[],
    uniform       real dataOut[])
{
  const uniform int ntasks = 128;
  launch[ntasks] laplace2d_task(
      ntasks,
      nx, ny, 
      dataIn, dataOut);
}
