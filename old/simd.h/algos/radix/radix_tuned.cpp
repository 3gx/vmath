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



/********** radix sort **********/

struct RadixSort
{
  private:
    enum 
    {
      NUMBITS = 8,
      NUMBUCKETS = (1<<NUMBITS),
    };

    int count;
    int blockDim;
    int gridDim;
    int numBlocks;

    unsigned int *sorted;
    int *excScanBlockPtr, *countsBlockPtr;

  public:

    int get_NUMBITS() const {return NUMBITS;}

    RadixSort(const int _count) : count(_count)
  {
#pragma omp parallel
#pragma omp master
    gridDim = omp_get_num_threads();

    if (0)
    {
      blockDim = std::max((count/gridDim/64) & -64, 64);  /* sandy bridge */
    }
    else
    {
      blockDim = std::max((count/gridDim/4) & -64, 64);   /* xeonphi */
    }
    numBlocks = gridDim;
    blockDim = (count+gridDim-1)/gridDim;

    const int ntmp = numBlocks * NUMBUCKETS;
    posix_memalign((void**)&sorted, 64, count*sizeof(int));
    posix_memalign((void**)&excScanBlockPtr, 64, ntmp*sizeof(int));
    posix_memalign((void**)& countsBlockPtr, 64, ntmp*sizeof(int));

    int (*excScanBlock)[NUMBUCKETS] = (int (*)[NUMBUCKETS])excScanBlockPtr;
    int (* countsBlock)[NUMBUCKETS] = (int (*)[NUMBUCKETS]) countsBlockPtr;

#pragma omp parallel
    {
      const int blockIdx = omp_get_thread_num();
      for(int block = blockIdx; block < numBlocks; block += gridDim)
#pragma simd
        for (int i = 0; i < NUMBUCKETS; i++)
          countsBlock[block][i] = excScanBlock[block][i] = 0;
#pragma omp for
      for (int i = 0; i < count; i++)
        sorted[i] = 0;
    }
  } 

    ~RadixSort()
    {
      free(sorted);
      free(excScanBlockPtr);
      free(countsBlockPtr);
    }

  private:

#define uniform
    void countPass(
        const uniform int keysAll[],
        const uniform int bit,
        const uniform int numElements,
        uniform int countsAll[],
        uniform int countsGlobal[])
    {
      const int blockIdx = omp_get_thread_num();
      const uniform int numBlocks = omp_get_num_threads();
#if 0
      const unsigned int *keys = keys1 + blockIdx*blockDim;
      const uniform int  blockDim = (count1 + numBlocks - 1) / numBlocks;
      int count = std::min(count1 - blockIdx*blockDim, blockDim);
      int *counts = counts1 + blockIdx*NUMBUCKETS;
#endif

  const uniform int  blockDim = (numElements + numBlocks - 1) / numBlocks;

  const uniform int mask = (1 << NUMBITS) - 1;

  const uniform int * uniform keys =   keysAll + blockIdx*blockDim;
  uniform int * uniform     counts = countsAll + blockIdx*NUMBUCKETS;
  const uniform int           count = std::min(numElements - blockIdx*blockDim, blockDim);

      // Compute the histogram of radix digits for this block only. This
      // corresponds exactly to the count kernel in the GPU implementation.
#pragma simd
      for (int i = 0; i < NUMBUCKETS; i++)
        counts[i] = 0;
      for(int i = 0; i < count; ++i) 
      {
        const int key = mask & ((unsigned int)keys[i] >> bit);
        counts[key]++;
      }
      for (int digit = 0; digit < NUMBUCKETS; digit++)
#pragma omp atomic
        countsGlobal[digit] += counts[digit];
    }

    void sortPass(
    uniform int keysAll[],
    uniform int sorted[],
    uniform int bit,
    uniform int numElements,
    uniform int digitOffsetsAll[],
    uniform int sharedCounts[])
    {
      const int blockIdx = omp_get_thread_num();
      const uniform int numBlocks = omp_get_num_threads();
  const uniform int  blockDim = (numElements + numBlocks - 1) / numBlocks;


  uniform int * uniform localCounts = sharedCounts + blockIdx*NUMBUCKETS;

  const uniform int keyIndex = blockIdx * blockDim;
  uniform int * uniform keys = keysAll + keyIndex;
  uniform int * uniform digitOffsets = digitOffsetsAll + blockIdx*NUMBUCKETS;
  const uniform int nloc = std::min(numElements - keyIndex, blockDim);


      // Compute the histogram of radix digits for this block only. This
      // corresponds exactly to the count kernel in the GPU implementation.
      const int NUMBUCKETS = 1<< NUMBITS;
      const unsigned int mask = NUMBUCKETS - 1;

      for(int i = 0; i < NUMBUCKETS; i++)
        localCounts[i] = 0;
#if 1
      for(int i = 0; i < nloc; i++)
      {
        // Extract the key 
        const int key = mask & (keys[i]>> bit);
        const int rel = localCounts[key];
        const int scatter = rel + digitOffsets[key];

        sorted[scatter] = keys[i];

        localCounts[key] = 1 + rel;
      }
#endif
    }

  public:
    void exclusive_scan(int *excScanPtr, int *countsBlockPtr, int numBlocks)
    {
      int (*excScanBlock)[NUMBUCKETS] = (int (*)[NUMBUCKETS])excScanBlockPtr;
      int (* countsBlock)[NUMBUCKETS] = (int (*)[NUMBUCKETS]) countsBlockPtr;
#if 1
#pragma omp parallel
      {
        const int blockIdx = omp_get_thread_num();
        for (int digit = blockIdx; digit < NUMBUCKETS; digit += gridDim)
        {
          int prev = excScanBlock[0][digit];
          for (int block = 0; block < numBlocks; block++)
          {
            const int y = countsBlock[block][digit];
            excScanBlock[block][digit] = prev;
            prev += y;
          }
        }
      }
#else
      int NTMAX = 1024;
      int partial_sum[NTMAX][NUMBUCKETS];
      assert(gridDim < NTMAX);
      for (int i = 0; i < NTMAX; i++)
        for (int j = 0; j < NUMBUCKETS; j++)
          partial_sum[i][j] = 0;

      const int blockDim = (numBlocks+gridDim-1)/gridDim;
      assert(blockDim == 1);
#pragma omp parallel
      {
        const int blockIdx = omp_get_thread_num();
        const int bbeg = blockIdx*blockDim;
        const int bend = std::min(bbeg+blockDim, numBlocks);

        for (int digit = 0; digit < NUMBUCKETS; digit++)
        {
          int prev = bbeg == 0 ? excScanBlock[0][digit] : 0;
          for (int block = bbeg; block < bend; block++)
          {
            const int y = countsBlock[block][digit];
            excScanBlock[block][digit] = prev;
            prev += y;
          }

          partial_sum[blockIdx][digit] = 
            excScanBlock[bend-1][digit] + countsBlock[bend-1][digit];
        }
      }

#pragma omp parallel
      {
        int prefix_sum[NTMAX];
        const int blockIdx = omp_get_thread_num();
        for (int digit = blockIdx; digit < NUMBUCKETS; digit += gridDim)
          if (digit < NUMBUCKETS) 
          {
            for (int i = 0; i < gridDim; i++)
              prefix_sum[i] = partial_sum[i][digit];

            int prev = 0;
            for (int i = 0; i < gridDim; i++)
            {
              const int y = prefix_sum[i];
              partial_sum[i][digit] = prev;
              prev += y;
            }
          }
      }

#pragma omp parallel
      {
        const int blockIdx = omp_get_thread_num();
        const int bbeg = blockIdx*blockDim;
        const int bend = std::min(bbeg+blockDim, numBlocks);
        for (int digit = 0; digit < NUMBUCKETS; digit++)
          for (int block = bbeg; block < bend; block++)
            excScanBlock[block][digit] += partial_sum[blockIdx][digit];
      }
#endif

    }

    void sort(unsigned int *keys)
    {
      int  countsGlobal[NUMBUCKETS] __attribute__((aligned(64))) = {0};

      int (*excScanBlock)[NUMBUCKETS] = (int (*)[NUMBUCKETS])excScanBlockPtr;
      int (* countsBlock)[NUMBUCKETS] = (int (*)[NUMBUCKETS]) countsBlockPtr;


#if 0
#define PROFILE
#endif

#ifdef PROFILE
      double dt1, dt2, dt3, dt4, dt5;
      double t0,t1;
      dt1=dt2=dt3=dt4=dt5=0.0;

      const double tbeg = rtc();
#endif


      for(int bit = 0; bit < 32; bit += NUMBITS)
      {

#ifdef PROFILE
        t0 = rtc();
#endif
        
        for (int digit = 0; digit < NUMBUCKETS; digit++)
          countsGlobal[digit] = 0;

#pragma omp parallel
        {
          /* histogramming each of the block */
            countPass(
                (int*)keys,
                bit,
                count,
                (int*)countsBlock,
                countsGlobal);
        }



#ifdef PROFILE
        { t1 = rtc(); dt1 += t1 - t0; t0 = t1; }
#endif



#ifdef PROFILE
        { t1 = rtc(); dt2 += t1 - t0; t0 = t1; }
#endif

        /* exclusive scan on the histogram */
        excScanBlock[0][0] = 0;
        for(int digit = 1; digit < NUMBUCKETS; digit++)
          excScanBlock[0][digit] = excScanBlock[0][digit-1] + countsGlobal[digit-1];


#ifdef PROFILE
        { t1 = rtc(); dt3 += t1 - t0; t0 = t1; }
#endif

        exclusive_scan(&excScanBlock[0][0], &countsBlock[0][0], numBlocks);

#ifdef PROFILE
        { t1 = rtc(); dt4 += t1 - t0; t0 = t1; }
#endif

        int counts[NUMBUCKETS*1024] ;
#pragma omp parallel
        {

          /* sorting */
            sortPass(
                (int*)keys,
                (int*)sorted,
                bit, 
                count, 
                (int*)excScanBlock,
                (int*)counts);
        }


#ifdef PROFILE
        { t1 = rtc(); dt5 += t1 - t0; t0 = t1; }
#endif

        std::swap(keys, sorted);
      }

      const double tend = rtc();

#ifdef PROFILE
      printf("dt1= %g \n", dt1);
      printf("dt2= %g \n", dt2);
      printf("dt3= %g \n", dt3);
      printf("dt4= %g \n", dt4);
      printf("dt5= %g \n", dt5);
      printf("dt = %g \n", tend-tbeg);
#endif
    }

};

int main(int argc, char * argv[])
{
  const int n = argc > 1 ? atoi(argv[1]) : 10000000;

  unsigned int *data_in, *data_out, *data_tmp;
  posix_memalign((void**)&data_in,  64, n*sizeof(int));
  posix_memalign((void**)&data_out, 64, n*sizeof(int));
  posix_memalign((void**)&data_tmp, 64, n*sizeof(int));

  srand48(rtc()*65536);

#pragma omp parallel for
  for (int i = 0; i < n; i++)
  {
    data_in[i]  = drand48() * (1<<30);
  }
  std::random_shuffle(data_in, data_in + n);
#pragma omp parallel for
  for (int i = 0; i < n; i++)
  {
    data_out[i] = data_in[i];
    data_tmp[i] = data_in[i];
  }

  int nthreads = 0;
#pragma omp parallel
#pragma omp master
  nthreads = omp_get_num_threads();
  printf("exclusive_scan:: nthreads= %d\n", nthreads);

#if 0
  const double t0 = rtc();
  for (int i = 0; i < 1; i++)
    radixSort(data_out,  n);
  const double dt = rtc() - t0;
#else
  RadixSort r(n);
  const double t0 = rtc();
  for (int i = 0; i < 1; i++)
  {
    r.sort(data_out);
    //   std::sort(&data_out[0],&data_out[n]);
  }
  const double dt = rtc() - t0;
#endif
  printf("scan done in %g sec, rate= %g  MEL/s  BW= %g GBs \n",
      dt, n/dt/1e6, sizeof(int)*n*(32/r.get_NUMBITS())*3/dt/1e9);


#if 1
  std::vector<int> gold(n);
  for (int i = 0; i < n; i++)
    gold[i] = data_tmp[i];
  std::sort(&gold[0],&gold[n]);

  for (int i = 0; i < n; i++)
    assert(data_out[i] == gold[i]);
#endif

  free(data_out);

  return 0;
}

