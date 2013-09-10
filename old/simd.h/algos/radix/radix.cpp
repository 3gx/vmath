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
      PAD = 1,
      numBits = 8,
      numBuckets = (1<<numBits),
      numBucketsPad = numBuckets * PAD
    };

    int count;
    int blockSize;
    int gridDim;
    int numBlocks;

    unsigned int *sorted;
    int *excScanBlockPtr, *countsBlockPtr;

  public:

    int get_numBits() const {return numBits;}

    RadixSort(const int _count) : count(_count)
  {
#pragma omp parallel
#pragma omp master
    gridDim = omp_get_num_threads();

    if (0)
    {
      blockSize = std::max((count/gridDim/64) & -64, 64);  /* sandy bridge */
    }
    else
    {
      blockSize = std::max((count/gridDim/4) & -64, 64);   /* xeonphi */
    }

    numBlocks  = (count + blockSize - 1) / blockSize;

    const int ntmp = numBlocks * numBucketsPad;
    posix_memalign((void**)&sorted, 64, count*sizeof(int));
    posix_memalign((void**)&excScanBlockPtr, 64, ntmp*sizeof(int));
    posix_memalign((void**)& countsBlockPtr, 64, ntmp*sizeof(int));

    int (*excScanBlock)[numBucketsPad] = (int (*)[numBucketsPad])excScanBlockPtr;
    int (* countsBlock)[numBucketsPad] = (int (*)[numBucketsPad]) countsBlockPtr;

#pragma omp parallel
    {
      const int blockIdx = omp_get_thread_num();
      for(int block = blockIdx; block < numBlocks; block += gridDim)
#pragma simd
        for (int i = 0; i < numBuckets; i++)
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

    void countPass(
        const unsigned int* keys, 
        const   int bit,
        const int count,
        int* counts) 
    {
      // Compute the histogram of radix digits for this block only. This
      // corresponds exactly to the count kernel in the GPU implementation.
      const unsigned int mask = (1 << numBits) - 1;
#pragma simd
      for (int i = 0; i < numBuckets; i++)
        counts[i] = 0;
#if 1
      for(int i = 0; i < count; ++i) 
      {
        const int key = mask & (keys[i] >> bit);
        counts[key]++;
      }
#endif
    }

    void sortPass(
        const unsigned int* keys,
        unsigned int* sorted, 
        int bit, 
        const int count,
        const int* digitOffsets,
        int* counts)
    {

      // Compute the histogram of radix digits for this block only. This
      // corresponds exactly to the count kernel in the GPU implementation.
      const int numBuckets = 1<< numBits;
      const unsigned int mask = numBuckets - 1;

#if 1
      for(int i = 0; i < count; i++)
      {
        // Extract the key 
        const int key = mask & (keys[i]>> bit);
        const int rel = counts[key];
        const int scatter = rel + digitOffsets[key];

        sorted[scatter] = keys[i];

        counts[key] = 1 + rel;
      }
#endif
    }

  public:

    void sort(unsigned int *keys)
    {
      int  countsGlobal[numBuckets] __attribute__((aligned(64))) = {0};
      int excScanGlobal[numBuckets] __attribute__((aligned(64))) = {0};
      int excScanBlockL[numBuckets] __attribute__((aligned(64))) = {0};
      int  digitOffsets[numBuckets] __attribute__((aligned(64))) = {0};

      int (*excScanBlock)[numBucketsPad] = (int (*)[numBucketsPad])excScanBlockPtr;
      int (* countsBlock)[numBucketsPad] = (int (*)[numBucketsPad]) countsBlockPtr;


#if 0
#define PROFILE
#endif

#ifdef PROFILE
      double dt1, dt2, dt3, dt4, dt5;
      double t0,t1;
      dt1=dt2=dt3=dt4=dt5=0.0;

      const double tbeg = rtc();
#endif

#pragma omp parallel
      {
        const int blockIdx = omp_get_thread_num();

        for(int bit = 0; bit < 32; bit += numBits)
        {

#ifdef PROFILE
#pragma omp master
          t0 = rtc();
#endif

          /* histogramming each of the block */
          for(int block = blockIdx; block < numBlocks; block += gridDim)
            countPass(
                keys + block*blockSize,
                bit,
                std::min(count - block*blockSize, blockSize),
                &countsBlock[block][0]);


#pragma omp barrier

#ifdef PROFILE
#pragma omp master
          { t1 = rtc(); dt1 += t1 - t0; t0 = t1; }
#endif


          /* compute global histogram */
          for (int digit = blockIdx; digit < numBuckets; digit += gridDim)
          {
            int sum = 0.0;
            for (int block = 0; block < numBlocks; block++)
              sum += countsBlock[block][digit];
            countsGlobal[digit] = sum;
          }

#pragma omp barrier

#ifdef PROFILE
#pragma omp master
          { t1 = rtc(); dt2 += t1 - t0; t0 = t1; }
#endif

          /* exclusive scan on the histogram */
#pragma omp single
            for(int digit = 1; digit < numBuckets; digit++)
                excScanGlobal[digit] = excScanGlobal[digit - 1] + countsGlobal[digit - 1];

#pragma omp barrier

#ifdef PROFILE
#pragma omp master
          { t1 = rtc(); dt3 += t1 - t0; t0 = t1; }
#endif

          /* computing offsets for each digit */
          for (int digit = blockIdx; digit < numBuckets; digit += gridDim)
          {
            int dgt =  digitOffsets[digit];
            for (int block = 0; block < numBlocks; block++)
            {
              excScanBlock[block][digit] = dgt + excScanGlobal[digit];
              dgt += countsBlock[block][digit];
            }
            digitOffsets[digit] = 0;
          }

#pragma omp barrier

#ifdef PROFILE
#pragma omp master
          { t1 = rtc(); dt4 += t1 - t0; t0 = t1; }
#endif

          /* sorting */
          for(int block = blockIdx; block < numBlocks; block += gridDim)
          {
            int counts[numBuckets] = {0};

            const int keyIndex = block * blockSize;
            sortPass(
                keys + keyIndex, 
                sorted,
                bit, 
                std::min(count - keyIndex, blockSize), 
                &excScanBlock[block][0],
                &counts[0]);
          }

#pragma omp barrier

#ifdef PROFILE
#pragma omp master
          { t1 = rtc(); dt5 += t1 - t0; t0 = t1; }
#endif

#pragma omp single
          {
#pragma simd
            for (int i = 0; i < numBuckets; i++)
              countsGlobal[i] = excScanGlobal[i]  = 0;
            std::swap(keys, sorted);
          }
        }
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
  const int n = argc > 1 ? atoi(argv[1]) : 112000000;


  unsigned int *data_in, *data_out, *data_tmp;
  posix_memalign((void**)&data_in,  64, n*sizeof(int));
  posix_memalign((void**)&data_out, 64, n*sizeof(int));
  posix_memalign((void**)&data_tmp, 64, n*sizeof(int));

  srand48(rtc()*65536);

#pragma omp parallel for
  for (int i = 0; i < n; i++)
  {
    data_in[i]  = drand48() * (1<<30);
    data_out[i] = data_in[i];
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
      dt, n/dt/1e6, sizeof(int)*n*(32/r.get_numBits())*3/dt/1e9);


#if 1
  std::vector<int> gold(n);
  for (int i = 0; i < n; i++)
    gold[i] = data_in[i];
  std::sort(&gold[0],&gold[n]);

  for (int i = 0; i < n; i++)
    assert(data_out[i] == gold[i]);
#endif

  free(data_out);

  return 0;
}

