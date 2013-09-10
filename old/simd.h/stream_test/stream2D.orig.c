#include <assert.h>
/*-----------------------------------------------------------------------*/
/* Program: Stream                                                       */
/* Revision: $Id: stream.c,v 5.9 2009/04/11 16:35:00 mccalpin Exp $ */
/* Original code developed by John D. McCalpin                           */
/* Programmers: John D. McCalpin                                         */
/*              Joe R. Zagar                                             */
/*                                                                       */
/* This program measures memory transfer rates in MB/s for simple        */
/* computational kernels coded in C.                                     */
/*-----------------------------------------------------------------------*/
/* Copyright 1991-2005: John D. McCalpin                                 */
/*-----------------------------------------------------------------------*/
/* License:                                                              */
/*  1. You are free to use this program and/or to redistribute           */
/*     this program.                                                     */
/*  2. You are free to modify this program for your own use,             */
/*     including commercial use, subject to the publication              */
/*     restrictions in item 3.                                           */
/*  3. You are free to publish results obtained from running this        */
/*     program, or from works that you derive from this program,         */
/*     with the following limitations:                                   */
/*     3a. In order to be referred to as "STREAM benchmark results",     */
/*         published results must be in conformance to the STREAM        */
/*         Run Rules, (briefly reviewed below) published at              */
/*         http://www.cs.virginia.edu/stream/ref.html                    */
/*         and incorporated herein by reference.                         */
/*         As the copyright holder, John McCalpin retains the            */
/*         right to determine conformity with the Run Rules.             */
/*     3b. Results based on modified source code or on runs not in       */
/*         accordance with the STREAM Run Rules must be clearly          */
/*         labelled whenever they are published.  Examples of            */
/*         proper labelling include:                                     */
/*         "tuned STREAM benchmark results"                              */
/*         "based on a variant of the STREAM benchmark code"             */
/*         Other comparable, clear and reasonable labelling is           */
/*         acceptable.                                                   */
/*     3c. Submission of results to the STREAM benchmark web site        */
/*         is encouraged, but not required.                              */
/*  4. Use of this program or creation of derived works based on this    */
/*     program constitutes acceptance of these licensing restrictions.   */
/*  5. Absolutely no warranty is expressed or implied.                   */
/*-----------------------------------------------------------------------*/
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

#ifndef N
#define NX 2048
#endif
#ifndef NTIMES
#   define NTIMES	10
#endif
#ifndef OFFSET
#   define OFFSET 8
#endif

/*
 *	3) Compile the code with full optimization.  Many compilers
 *	   generate unreasonably bad code before the optimizer tightens
 *	   things up.  If the results are unreasonably good, on the
 *	   other hand, the optimizer might be too smart for me!
 *
 *         Try compiling with:
 *               cc -O stream_omp.c -o stream_omp
 *
 *         This is known to work on Cray, SGI, IBM, and Sun machines.
 *
 *
 *	4) Mail the results to mccalpin@cs.virginia.edu
 *	   Be sure to include:
 *		a) computer hardware model number and software revision
 *		b) the compiler flags
 *		c) all of the output from the test case.
 * Thanks!
 *
 */

# define HLINE "-------------------------------------------------------------\n"

# ifndef MIN
# define MIN(x,y) ((x)<(y)?(x):(y))
# endif
# ifndef MAX
# define MAX(x,y) ((x)>(y)?(x):(y))
# endif

static double 
aa[NX][NX+OFFSET],
  bb[NX][NX+OFFSET], 
  cc[NX][NX+OFFSET],
  dd[NX][NX+OFFSET], 
  ee[NX][NX+OFFSET], 
  ff[NX][NX+OFFSET],
  gg[NX][NX+OFFSET], 
  hh[NX][NX+OFFSET];

typedef struct 
{
  double  cc;
#if 1
  double aa;
  double bb;
#endif
#if 1
  double dd;
  double ee;
#endif
  double  ff;
#if 1
  double gg;
  double hh;
#endif
} mydataS;

static mydataS mydata[NX][NX+OFFSET];

static double	avgtime[4] = {0}, maxtime[4] = {0},
							mintime[4] = {FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};

static char	*label[4] = {"Copy:      ", "Scale:     ",
	"Add:       ", "Triad:     "};

static double	bytes[4] = {
	8 * sizeof(double) * NX*NX,
	2 * sizeof(double) * NX*NX,
	3 * sizeof(double) * NX*NX,
	3 * sizeof(double) * NX*NX,
};

extern double mysecond();
extern void checkSTREAMresults();
#ifdef TUNED
extern void tuned_STREAM_Copy();
extern void tuned_STREAM_Scale(double scalar);
extern void tuned_STREAM_Add();
extern void tuned_STREAM_Triad(double scalar);
#endif
#ifdef _OPENMP
extern int omp_get_num_threads();
#endif
	int
main()
{
	int			quantum, checktick();
	int			BytesPerWord;
	register int	i, j, k;
	double		scalar, t, times[4][NTIMES];

	/* --- SETUP --- determine precision and check timing --- */

	printf(HLINE);
	printf("STREAM version $Revision: 5.9 $\n");
	printf(HLINE);
	BytesPerWord = sizeof(double);
	printf("This system uses %d bytes per DOUBLE PRECISION word.\n",
			BytesPerWord);

	printf(HLINE);
#ifdef NO_LONG_LONG
	printf("Array size = %d, Offset = %d\n" , NX*NX, OFFSET);
#else
	printf("Array size = %llu, Offset = %d\n", (unsigned long long) NX*NX, OFFSET);
#endif

	printf("Total memory required = %.1f MB.\n",
			(3.0 * BytesPerWord) * ( (double) NX*NX / 1048576.0));
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

	/* Get initial value for system clock. */
#if 0
#pragma omp parallel for
	for (j=0; j<NX; j++)
    for (i = 0; i < NX; i++)
    {
      aa[j][i] = 1.0;
      bb[j][i] = 2.0;
      cc[j][i] = 0.0;
    }
#endif

  printf(HLINE);

  if  ( (quantum = checktick()) >= 1) 
    printf("Your clock granularity/precision appears to be "
        "%d microseconds.\n", quantum);
  else {
    printf("Your clock granularity appears to be "
        "less than one microsecond.\n");
    quantum = 1;
  }

  t = mysecond();
#if 0
#pragma omp parallel for
  for (j = 0; j < NX; j++)
    for (i = 0; i < NX; i++)
    aa[j][i] = 2.0E0 * aa[j][i];
#endif
  t = 1.0E6 * (mysecond() - t);

  printf("Each test below will take on the order"
      " of %d microseconds.\n", (int) t  );
  printf("   (= %d clock ticks)\n", (int) (t/quantum) );
  printf("Increase the size of the arrays if this shows that\n");
  printf("you are not getting at least 20 clock ticks per test.\n");

  printf(HLINE);

  printf("WARNING -- The above is only a rough guideline.\n");
  printf("For best results, please be sure you know the\n");
  printf("precision of your system timer.\n");
  printf(HLINE);

  /*	--- MAIN LOOP --- repeat test cases NTIMES times --- */


  scalar = 3.0;
  for (k=0; k<NTIMES; k++)
  {
    times[0][k] = mysecond();
#ifdef TUNED
    tuned_STREAM_Copy();
#else
#pragma omp parallel for
    for (j=0; j<NX; j++)
      for (i=0;i<NX;i++)
      cc[j][i] = aa[j][i];
#endif
    times[0][k] = mysecond() - times[0][k];

#if 0
    times[1][k] = mysecond();
#ifdef TUNED
    tuned_STREAM_Scale(scalar);
#else
#pragma omp parallel for
    for (j=0; j<NX; j++)
      for (i=0;i<NX;i++)
      bb[j][i] = scalar*cc[j][i];
#endif
    times[1][k] = mysecond() - times[1][k];
    times[2][k] = mysecond();
#ifdef TUNED
    tuned_STREAM_Add();
#else
#pragma omp parallel for
    for (j=0; j<NX; j++)
      for (i=0;i<NX;i++)
      cc[j][i] = aa[j][i]+bb[j][i];
#endif
    times[2][k] = mysecond() - times[2][k];

    times[3][k] = mysecond();
#ifdef TUNED
    tuned_STREAM_Triad(scalar);
#else
#pragma omp parallel for
    for (j=0; j<NX; j++)
      for (i=0;i<NX;i++)
      aa[j][i] = bb[j][i]+scalar*cc[j][i];
#endif
    times[3][k] = mysecond() - times[3][k];
#endif
  }

  /*	--- SUMMARY --- */

  for (k=1; k<NTIMES; k++) /* note -- skip first iteration */
  {
    for (j=0; j<4; j++)
    {
      avgtime[j] = avgtime[j] + times[j][k];
      mintime[j] = MIN(mintime[j], times[j][k]);
      maxtime[j] = MAX(maxtime[j], times[j][k]);
    }
  }

  printf("Function      Rate (MB/s)   Avg time     Min time     Max time\n");
  for (j=0; j<4; j++) {
    avgtime[j] = avgtime[j]/(double)(NTIMES-1);

    printf("%s%11.4f  %11.4f  %11.4f  %11.4f\n", label[j],
        1.0E-06 * bytes[j]/mintime[j],
        avgtime[j],
        mintime[j],
        maxtime[j]);
  }
  printf(HLINE);

  /* --- Check Results --- */
  checkSTREAMresults();
  printf(HLINE);

  return 0;
}

# define	M	20

  int
checktick()
{
  int		i, minDelta, Delta;
  double	t1, t2, timesfound[M];

  /*  Collect a sequence of M unique time values from the system. */

  for (i = 0; i < M; i++) {
    t1 = mysecond();
    while( ((t2=mysecond()) - t1) < 1.0E-6 )
      ;
    timesfound[i] = t1 = t2;
  }

  /*
   * Determine the minimum difference between these M values.
   * This result will be our estimate (in microseconds) for the
   * clock granularity.
   */

  minDelta = 1000000;
  for (i = 1; i < M; i++) {
    Delta = (int)( 1.0E6 * (timesfound[i]-timesfound[i-1]));
    minDelta = MIN(minDelta, MAX(Delta,0));
  }

  return(minDelta);
}



/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.  */

#include <sys/time.h>

double mysecond()
{
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

void checkSTREAMresults ()
{
  double aj,bj,cj,scalar;
  double asum,bsum,csum;
  double epsilon;
  int	i,j,k;

  /* reproduce initialization */
  aj = 1.0;
  bj = 2.0;
  cj = 0.0;
  /* a[] is modified during timing check */
  aj = 2.0E0 * aj;
  /* now execute timing loop */
  scalar = 3.0;
  for (k=0; k<NTIMES; k++)
  {
    cj = aj;
    bj = scalar*cj;
    cj = aj+bj;
    aj = bj+scalar*cj;
  }
  aj = aj * (double) (NX*NX);
  bj = bj * (double) (NX*NX);
  cj = cj * (double) (NX*NX);

  asum = 0.0;
  bsum = 0.0;
  csum = 0.0;
  for (j=0; j<NX; j++) 
  for (i=0; i<NX; i++) {
    asum += aa[j][i];
    bsum += bb[j][i];
    csum += cc[j][i];
  }
#ifdef VERBOSE
  printf ("Results Comparison: \n");
  printf ("        Expected  : %f %f %f \n",aj,bj,cj);
  printf ("        Observed  : %f %f %f \n",asum,bsum,csum);
#endif

#ifndef abs
#define abs(a) ((a) >= 0 ? (a) : -(a))
#endif
  epsilon = 1.e-8;

  if (abs(aj-asum)/asum > epsilon) {
    printf ("Failed Validation on array a[]\n");
    printf ("        Expected  : %f \n",aj);
    printf ("        Observed  : %f \n",asum);
  }
  else if (abs(bj-bsum)/bsum > epsilon) {
    printf ("Failed Validation on array b[]\n");
    printf ("        Expected  : %f \n",bj);
    printf ("        Observed  : %f \n",bsum);
  }
  else if (abs(cj-csum)/csum > epsilon) {
    printf ("Failed Validation on array c[]\n");
    printf ("        Expected  : %f \n",cj);
    printf ("        Observed  : %f \n",csum);
  }
  else {
    printf ("Solution Validates\n");
  }
}

void tuned_STREAM_Copy()
{
  int i,j;
#if 1
#pragma omp parallel for
  for (j=1; j<NX-1; j++)
    for (i=1; i<NX-1; i++)
    {
#if 1
      cc[j][i] = 0
#if 1
        + aa[j][i-1]+aa[j  ][i+1] + aa[j][i]  
        + aa[j-1][i]+aa[j-1][i  ] + aa[j][i] 
#endif
        + bb[j][i-1]+bb[j  ][i+1] + bb[j][i] 
        + bb[j-1][i]+bb[j-1][i  ] + bb[j][i]
        + dd[j][i-1]+dd[j  ][i+1] + dd[j][i]  
        + dd[j-1][i]+dd[j-1][i  ] + dd[j][i] 
        + ee[j][i-1]+ee[j  ][i+1] + ee[j][i]  
        + ee[j-1][i]+ee[j-1][i  ] + ee[j][i] 
        + ff[j][i-1]+ff[j  ][i+1] + ff[j][i]  
        + ff[j-1][i]+ff[j-1][i  ] + ff[j][i]
#if 1
        + gg[j][i-1]+gg[j  ][i+1] + gg[j][i]  
        + gg[j-1][i]+gg[j-1][i  ] + gg[j][i] 
        + hh[j][i-1]+hh[j  ][i+1] + hh[j][i]  
        + hh[j-1][i]+hh[j-1][i  ] + hh[j][i]
#endif
        ;
#else
      mydata[j][i].cc = 0
#if 1
        + mydata[j][i-1].aa+mydata[j  ][i+1].aa + mydata[j][i].aa
        + mydata[j-1][i].aa+mydata[j-1][i  ].aa + mydata[j][i].aa
        + mydata[j][i-1].bb+mydata[j  ][i+1].bb + mydata[j][i].bb 
        + mydata[j-1][i].bb+mydata[j-1][i  ].bb + mydata[j][i].bb
        + mydata[j][i-1].dd+mydata[j  ][i+1].dd + mydata[j][i].dd  
        + mydata[j-1][i].dd+mydata[j-1][i  ].dd + mydata[j][i].dd 
        + mydata[j][i-1].ee+mydata[j  ][i+1].ee + mydata[j][i].ee  
        + mydata[j-1][i].ee+mydata[j-1][i  ].ee + mydata[j][i].ee 
#endif
        + mydata[j][i-1].ff+mydata[j  ][i+1].ff + mydata[j][i].ff  
        + mydata[j-1][i].ff+mydata[j-1][i  ].ff + mydata[j][i].ff
#if 1
        + mydata[j][i-1].gg+mydata[j  ][i+1].gg + mydata[j][i].gg
        + mydata[j-1][i].gg+mydata[j-1][i  ].gg + mydata[j][i].gg
        + mydata[j][i-1].hh+mydata[j  ][i+1].hh + mydata[j][i].hh 
        + mydata[j-1][i].hh+mydata[j-1][i  ].hh + mydata[j][i].hh
#endif
#endif
        ;
    }
#else
#pragma omp parallel
  {
#pragma omp for
    for (j=1; j<NX-1; j++)
      for (i=1; i<NX-1; i++)
      {
        cc[j][i] = 0
          + aa[j][i-1]+aa[j  ][i+1] + aa[j][i]  
          + aa[j-1][i]+aa[j-1][i  ] + aa[j][i] 
          + bb[j][i-1]+bb[j  ][i+1] + bb[j][i] 
          + bb[j-1][i]+bb[j-1][i  ] + bb[j][i]
          + dd[j][i-1]+dd[j  ][i+1] + dd[j][i]  
          + dd[j-1][i]+dd[j-1][i  ] + dd[j][i];
      }

#pragma omp for
    for (j=1; j<NX-1; j++)
      for (i=1; i<NX-1; i++)
      {
        cc[j][i] +=
          + ee[j][i-1]+ee[j  ][i+1] + ee[j][i]  
          + ee[j-1][i]+ee[j-1][i  ] + ee[j][i] 
          + ff[j][i-1]+ff[j  ][i+1] + ff[j][i]  
          + ff[j-1][i]+ff[j-1][i  ] + ff[j][i]
          + gg[j][i-1]+gg[j  ][i+1] + gg[j][i]  
          + gg[j-1][i]+gg[j-1][i  ] + gg[j][i] 
          + hh[j][i-1]+hh[j  ][i+1] + hh[j][i]  
          + hh[j-1][i]+hh[j-1][i  ] + hh[j][i];
      }
  }
#endif

}

void tuned_STREAM_Scale(double scalar)
{
  int i,j;
#pragma omp parallel for
  for (j=0; j<NX; j++)
    for (i=0; i<NX; i++)
      bb[j][i] = scalar*cc[j][i];
}

void tuned_STREAM_Add()
{
  int i,j;
#pragma omp parallel for
  for (j=0; j<NX; j++)
    for (i=0; i<NX; i++)
    {
      cc[j][i] = aa[j][i]+bb[j][i];
    }
}

void tuned_STREAM_Triad(double scalar)
{
  int i, j;
#pragma omp parallel for
  for (j=0; j<NX; j++)
    for (i=0; i<NX; i++)
      aa[j][i] = bb[j][i]+scalar*cc[j][i];
}
