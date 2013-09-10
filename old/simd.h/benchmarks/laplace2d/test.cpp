#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "mytimer.h";

int main(int argc, char * argv[])
{
  assert(argc > 2);
  const int  n = atoi(argv[1]);
  fprintf(stderr, " n = %d\n", n);
  const double  c = atoi(argv[2]);
  fprintf(stderr, " c = %g\n", c);

  double *data = (double*)malloc(sizeof(double)*n);

#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; i++)
    data[i] = drand48();

  const double t0 = get_wtime();

  const int nrep = 10;
#pragma omp parallel for schedule(static)
  for (int k = 0; k < nrep; k++)
    for (int i = 1 ; i < n; i++)
      data[i] = data[i-1] * c;


  const double t1 = get_wtime();

  double res = 0.0;
  for (int i = 0; i < n; i++)
    res += data[i];

  fprintf(stderr, " res= %g  done in %g sec\n", res, t1-t0);


  return 0;
};
