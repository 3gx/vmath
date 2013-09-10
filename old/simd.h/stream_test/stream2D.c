#include <assert.h>

# include <stdio.h>
# include <math.h>
# include <float.h>
# include <limits.h>
# include <sys/time.h>
#include <math.h>
#include <stdlib.h>

/* INSTRUCTIONS:
 *
 *	1) Stream requires a good bit of memory to run.  Adjust the
 *          value of 'N' (below) to give a 'timing calibration' of 
 *          at least 20 clock-ticks.  This will provide rate estimates
 *          that should be good to about 5% precision.
 */

#ifndef NX
#define NX 1024
#endif
#ifndef ND
#define ND 8
#endif
#ifndef NTIMES
#define NTIMES 5
#endif
#ifndef OFFSET
#define OFFSET 8
#endif

# define HLINE "-------------------------------------------------------------\n"

# ifndef MIN
# define MIN(x,y) ((x)<(y)?(x):(y))
# endif
# ifndef MAX
# define MAX(x,y) ((x)>(y)?(x):(y))
# endif

static char	*label[2] = {"Array:      ", "Struct:     "};
static double	bytes[2] = {
	ND * sizeof(double) * NX*(NX+OFFSET),
	ND * sizeof(double) * NX*(NX+OFFSET)
};

double mysecond();


static double data[ND][NX][NX+OFFSET];
static struct
{
  double data[ND];
} list[NX][NX+OFFSET];

# define HLINE "-------------------------------------------------------------\n"

static double	avgtime[2] = {0}, maxtime[2] = {0},
							mintime[2] = {FLT_MAX,FLT_MAX};

int main()
{
	int			BytesPerWord;
	register int	i, j, k, kk;
	double		scalar, t, times[2][NTIMES];

  printf(HLINE);
  printf("ND= %d\n", ND);
  printf("NX= %d\n", NX);
  printf(HLINE);
  

	BytesPerWord = sizeof(double);
	printf("This system uses %d bytes per DOUBLE PRECISION word.\n",
			BytesPerWord);

	printf("Total memory required = %.1f MB.\n",
			(ND * BytesPerWord) * ( (double) NX*(NX+OFFSET) / 1048576.0));
	printf("Each test is run %d times, but only\n", NTIMES);
	printf("the *best* time for each is used.\n");

#ifdef _OPENMP
	printf(HLINE);
#pragma omp parallel 
	{
#pragma omp master
		{
			k = omp_get_num_threads();
			printf ("Number of Threads requested = %i\n",k);
		}
	}
#endif

	printf(HLINE);
#pragma omp parallel
	{
		printf ("Printing one line per active thread....\n");
	}

  printf(HLINE);

  /*	--- MAIN LOOP --- repeat test cases NTIMES times --- */


  for (k=0; k<NTIMES; k++)
  {
    times[0][k] = mysecond();
#pragma omp parallel for
    for (j = 1; j < NX-1; j++)
     for (i = 1; i < NX-1; i++)
     {
       double c;
       for (kk = 1; kk < ND; kk++)
        c += data[kk][j][i] +
          data[kk][j][i-1]+data[kk][j  ][i+1] + 
          data[kk][j-1][i]+data[kk][j-1][i  ];
       data[0][j][i] = c;
     }
    times[0][k] = mysecond() - times[0][k];

    times[1][k] = mysecond();
#pragma omp parallel for
    for (j = 1; j < NX-1; j++)
     for (i = 1; i < NX-1; i++)
     {
       double c;
       for (kk = 1; kk < ND; kk++)
        c += list[j][i].data[kk] +
          list[j][i-1].data[kk]+list[j  ][i+1].data[kk] + 
          list[j-1][i].data[kk]+list[j-1][i  ].data[kk];
       list[j][i].data[0] = c;
     }
    times[1][k] = mysecond() - times[1][k];
  }

  /*	--- SUMMARY --- */

  for (k=1; k<NTIMES; k++) /* note -- skip first iteration */
  {
		for (j=0; j<2; j++)
		{
			avgtime[j] =     avgtime[j] + times[j][k];
			mintime[j] = MIN(mintime[j],  times[j][k]);
			maxtime[j] = MAX(maxtime[j],  times[j][k]);
		}
  }

  printf("Function      Rate (MB/s)   Avg time     Min time     Max time\n");
  for (j=0; j<2; j++) {
    avgtime[j] = avgtime[j]/(double)(NTIMES-1);

    printf("%s%11.4f  %11.4f  %11.4f  %11.4f\n", label[j],
        1.0E-06 * bytes[j]/mintime[j],
        avgtime[j],
        mintime[j],
        maxtime[j]);
  }
  printf(HLINE);

  /* --- Check Results --- */
  printf(HLINE);

  return 0;
}


#include <sys/time.h>

double mysecond()
{
	struct timeval tp;
	struct timezone tzp;
	int i;

	i = gettimeofday(&tp,&tzp);
	return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
