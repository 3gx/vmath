/* segmented scan */
//    partial_sum[tx] = data[iend-1];


#include<cstdio>
#include<cstdlib>
#include<cassert>
#include<vector>
#include<omp.h>
#include<algorithm>

#include <sys/time.h>
inline double rtc() {
  struct timeval Tvalue;
  struct timezone dummy;

  gettimeofday(&Tvalue,&dummy);
  return ((double) Tvalue.tv_sec +1.e-6*((double) Tvalue.tv_usec));
}

#include <exclusive_scan_ispc.h>

static void exclusive_scan(const int n, int *data)
{
  const int NTHREAD_MAX = 65536;
  static int partial_sum[NTHREAD_MAX];

  static int nthreads = 0;
#pragma omp parallel
#pragma omp master
  nthreads = omp_get_num_threads();
  printf("exclusive_scan:: nthreads= %d\n", nthreads);

  const int nlaunch = std::min(n/16,nthreads);

#pragma omp parallel num_threads(nlaunch)
  {
    const int blockIdx = omp_get_thread_num();
    const int gridDim  = omp_get_num_threads();
    const int blockDim = std::max((n/gridDim) & -64, 16);


    if (blockIdx == 0)
    assert(gridDim < NTHREAD_MAX);

    int nblock = 0;
    for (int ibeg = blockIdx*blockDim; ibeg < n; ibeg += blockDim*gridDim, nblock += gridDim)
    {
      assert(nblock < NTHREAD_MAX);
      const int iend = std::min(ibeg+blockDim, n);

      const int value = data[iend-1];
#if 0
      ispc::exclusive_scan(iend-ibeg, &data[ibeg]);
#else
      int prev = 0;
      for (int i = ibeg; i < iend; i++)
      {
        const int y = data[i];
        data[i] = prev;
        prev += y;
      }
#endif
      partial_sum[nblock + blockIdx] = value + data[iend-1];
    }

#pragma omp barrier


#if 0
    if (blockIdx == 0)
        ispc::exclusive_scan(nblock, &partial_sum[0]);
#else
    if (blockIdx == 0)
    {
      int prev = 0;
      for (int i = 0; i < nblock; i++)
      {
        const int y = partial_sum[i];
        partial_sum[i] = prev;
        prev += y;
      }
    }
#endif

#pragma omp barrier

    nblock = 0;
    for (int ibeg = blockIdx*blockDim; ibeg < n; ibeg += blockDim*gridDim, nblock += gridDim)
    {
      const int iend = std::min(ibeg+blockDim, n);
#if 0
      ispc::add(iend-ibeg, partial_sum[nblock + blockIdx], &data[ibeg]);
#else  /* this one is slower */
      const int sum = partial_sum[nblock + blockIdx];
      for (int i = ibeg; i < iend; i++)
        data[i] += sum;
#endif
    }
  }
}

int main(int argc, char * argv[])
{
  const int n = argc > 1 ? atoi(argv[1]) : 112000000;


  int *data_in, *data_out;
  posix_memalign((void**)&data_in,  64, n*sizeof(int));
  posix_memalign((void**)&data_out, 64, n*sizeof(int));


#pragma omp parallel for
  for (int i = 0; i < n; i++)
  {
    data_in[i]  = i;
    data_out[i] = data_in[i];
  }

  const double t0 = rtc();
  exclusive_scan(n, &data_out[0]);
  const double dt = rtc() - t0;
  printf("scan done in %g sec, rate= %g  MEL/s  BW= %g GBs \n",
      dt, n/dt/1e6, sizeof(int)*n*4/dt/1e9);


  std::vector<int> gold(n);
  gold[0] = 0;
  for (int i = 0; i < n-1; i++)
    gold[i+1] = gold[i] + data_in[i];

  for (int i = 0; i < n; i++)
  {
//        printf("i= %d: orig= %d data= %d  gold= %d\n", i, data_in[i], data_out[i], gold[i]);
    assert(data_out[i] == gold[i]);

  }

  free(data_out);

  return 0;
}

