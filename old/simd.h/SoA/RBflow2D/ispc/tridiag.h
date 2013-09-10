#pragma once

template<typename VEC>
struct Tridiag
{
  VEC a, b, c, f;
  Tridiag() {}
  Tridiag(const VEC _a, const VEC _b, const VEC _c) :
    a(_a), b(_b), c(_c) {}
 
 template<const int N> 
   inline static void ReduceN(const int n, Tridiag mat[])
   {
     for (int j = 0; j < N; j++)
       mat[j].b = VEC(1.0)/mat[j].b;

     for (int i = N; i < n*N; i += N)
       for (int j = 0; j < N; j++)
       {
         mat[i+j  ].f  = VEC(-1.0) *mat[i+j  ].a*mat[i+j-N].b;
         mat[i+j  ].b += mat[i+j].f*mat[i+j-N].c;
         mat[i+j  ].b  = VEC(1.0)/  mat[i+j  ].b;
         mat[i+j-N].c  =           -mat[i+j-N].c;
       }
   }

 template<const int N>
   inline static void SolveN(const int n, const Tridiag mat[], VEC d[])
   {
     for (int i = N; i < n*N; i += N)
       for (int j = 0; j < N; j++)
         d[i+j] += mat[i+j].f*d[i+j-N];
     for (int j = 0; j < N; j++)
       d[n*N-N+j] *= mat[n*N-N+j].b;
     for (int i = n*N-2*N; i >= 0; i -= N)
       for (int j = 0; j < N; j++)
         d[i+j] = (d[i+j] + mat[i+j].c*d[i+j+N])*mat[i+j].b;
   }

 inline static void Reduce(const int n, Tridiag mat[])
 {
   assert(mat[0].b != VEC(0.0));
   mat[0].b = VEC(1.0)/mat[0].b;
   for (int i = 1; i < n; i++)
   {
     mat[i  ].f  = VEC(-1.0)*mat[i].a*mat[i-1].b;
     mat[i  ].b += mat[i].f*mat[i-1].c;
#if 1
     assert(mat[i].b != VEC(0.0));
#endif
     mat[i  ].b  = VEC(1.0)/mat[i].b;
     mat[i-1].c  = -mat[i-1].c;
   }
 }

 inline static void Solve(const int n, const Tridiag mat[], VEC d[])
 {
   for (int i = 1; i < n; i++)
     d[i] += mat[ i ].f*d[i-1];
   d[n-1] *= mat[n-1].b;
   for (int i = n-2; i >= 0; i--)
     d[i] = (d[i] + mat[i].c*d[i+1])*mat[i].b;
 }


 template<const int N>
   inline static void Solve(const int n, const Tridiag<VEC> mat[], VEC d[]) 
   {
     int iN = N;
     for (int i = 1; i < n; i++, iN += N)
       for (int j = 0; j < N; j++)
         d[iN+j] += mat[ i ].f*d[iN-N+j];
     iN -= N;
     for (int j = 0; j < N; j++)
       d[iN+j] *= mat[n-1].b;
     iN -= N;
     for (int i = n-2; i >= 0; i--, iN -= N)
       for (int j = 0; j < N; j++)
         d[iN+j] = (d[iN+j] + mat[i].c*d[iN+N+j])*mat[i].b;
   }

 inline static void Reduce(const int n, Tridiag mat[], VEC u[])
 {
   const VEC gamma    = -mat[0].b;
   const VEC invGamma = VEC(1.0)/gamma;

   mat[ 0 ].b -= gamma;
   mat[n-1].b -= mat[n-1].c*mat[0].a * invGamma;

   assert(mat[0].b != VEC(0.0));

   mat[0].b = VEC(1.0)/mat[0].b;
   u [0]   = gamma;
   for (int i = 1; i < n; i++)
   {
     mat[i  ].f  = VEC(-1.0)*mat[i].a*mat[i-1].b;
     mat[i  ].b += mat[i].f*mat[i-1].c;
#if 1
     assert(mat[i].b != VEC(0.0));
#endif
     mat[i  ].b = VEC(1.0)/mat[i].b;
     mat[i-1].c = -mat[i-1].c;
     u  [i  ]   =  mat[i  ].f*u[i-1];
   }
   u[n-1] += mat[n-1].c;
   u[n-1] *= mat[n-1].b;
   for (int i = n-2; i >= 0; i--)
     u[i] = (u[i] + mat[i].c*u[i+1])*mat[i].b;
 }

 inline static void Solve(const int n, const Tridiag mat[], const VEC u[], VEC d[])
 {
   for (int i = 1; i < n; i++)
     d[i] += mat[ i ].f*d[i-1];
   d[n-1] *= mat[n-1].b;
   for (int i = n-2; i >= 0; i--)
     d[i] = (d[i] + mat[i].c*d[i+1])*mat[i].b;

   const VEC invGamma = VEC(-2.0)*mat[0].b;
   const VEC fact = VEC(-1.0) *
     (d[0] + d[n-1]*mat[0].a*invGamma) /
     (u[0] + u[n-1]*mat[0].a*invGamma + VEC(1.0));

   for (int i = 0; i < n; i++)
     d[i] += fact*u[i];
 }

 template<const int N>
   inline static void Solve(const int n, const Tridiag mat[], const VEC u[], VEC d[])
   {
     int iN = N;
     for (int i = 1; i < n; i++, iN += N)
       for (int j = 0; j < N; j++)
         d[iN+j] += mat[i].f*d[iN-N+j];
     iN -= N;
     for (int j = 0; j < N; j++)
       d[iN+j] *= mat[n-1].b;
     iN -= N;
     for (int i = n-2; i >= 0; i--, iN -= N)
       for (int j = 0; j < N; j++)
         d[iN+j] = (d[iN+j] + mat[i].c*d[iN+N+j])*mat[i].b;

     const VEC invGamma = VEC(-2.0)* mat[0].b;
     const VEC invC     = invGamma * mat[0].a;

     VEC fact[N];
     iN = n*N-N;
     for (int j = 0; j < N; j++)
       fact[j] = VEC(-1.0) *
         (d[0+j] + d[iN+j]*invC) /
         (u[0  ] + u[ n-1]*invC + VEC(1.0));

     iN = 0;
     for (int i = 0; i < n; i++, iN += N)
       for (int j = 0; j < N; j++)
         d[iN+j] += fact[j]*u[i];
   }


};


void tridiag(const int n,
    real a[], real b[], real c[], real d[], real u[], real gam[])
{
  assert(b[0] != 0.0);  /* solution does not exist */

  real ibet = 1.0/b[0];
  u[0] = d[0]*ibet;

  for (int j = 1; j < n; j++)
  {
    gam[j] = c[j-1]*ibet;
    const real bet = b[j] - a[j]*gam[j];
    assert(bet != 0.0);   /* solution does not exist */
    ibet = 1.0/bet;
    u[j] = (d[j] - a[j]*u[j-1])*ibet;
  }

  for (int j = n-2; j >= 0; j--)
    u[j] -= gam[j+1]*u[j+1];
}


#if 0
void tripvmy_line(const int n,
    real ami[], real aci[], real api[], real rrr[],  
    real   q[], real   s[], real fei[])
{
  const int n1i = 0;
  const int n1f = n;

  const int ia = n1i + 1;
  const int ii = n1i + n1f;

  /*   COEFFICIENTS FOR TRIDIAGONAL INVERSION   */

  /*  THE INVERSION STARTS   */

  q[n1i] = -api[n1i]/aci[n1i];
  s[n1i] = -ami[n1i]/aci[n1i];

  const real fn = rrr[n1f-1];
  rrr[n1i] = rrr[n1i]/aci[n1i];


  /*     forward elimination sweep     */
  for (int i = ia; i < n1f; i++)
  {
    const real p = 1.0/(aci[i] + ami[i]*q[i-1]);
    q  [i] = - api[i]*p;
    s  [i] = - ami[i]*s[i-1]*p;
    rrr[i] = (rrr[i] - ami[i]*rrr[i-1])*p;
  }

  /*     backward pass     */
  s  [n1f-1] = 1.0;
  fei[n1f-1] = 0.0;

  for (int l = ia; l < n1f; l++)
  {
    const int i = ii - l - 1;
    s  [i] +=   q[i] * s[i+1];
    fei[i]  = rrr[i] + q[i]*fei[i+1];
  }

  const int i = ii - n1f;
  rrr[n1f-1] = 
    (fn-api[i]*fei[n1i] - ami[i]*fei[n1f-2]) / 
    (   api[i]*  s[n1i] + ami[i]*  s[n1f-2] + aci[i]);

  /*     backward elimination pass      */
  for (int l = ia; l < n1f; l++)
  {
    const int i = ii - l - 1;
    rrr[i] = rrr[n1f-1]*s[i] + fei[i];
  }
}

#else

void tripvmy_line(const int n,
    real ami[], real aci[], real api[], real rrr[],  
    real   u[], real   s[], real fei[])
{

  /*     Solves for a vector rrr(1:n1f) the “cyclic” set of linear equations */

  //     Choose a pivot as -b(1)
  const real gamma = -aci[0];

  //     Set up the diagonal of the modified tridiagonal system.

  const real aci0 = aci[0];
  const real acin = aci[n-1];

  aci[ 0 ] -= gamma;
  aci[n-1] -= api[n-1]*ami[0]/gamma;

  /*  Set up the vector u  */

#if 0
  /*      Solve A · x =r  */
  tridiag(n, ami, aci, api, rrr, s, fei);
  for (int i = 0; i < n; i++)
    rrr[i] = s[i];

  u[ 0 ] = gamma;
  u[n-1] = api[n-1];
  for (int i = 0; i < n-1; i++)
    u[i] = 0.0;
  /*      Solve A · z = u. */
  tridiag(n, ami, aci, api, u, s, fei);
  for (int i = 0; i < n; i++)
    u[i] = s[i];
#else
  /* reduce to upper tridiagonal */
  assert(aci[0] != 0.0);
  fei[0] = 1.0/aci[0];
  u [0] = gamma;
  for (int i = 1; i < n; i++)
  {
    assert(aci[i] != 0.0);
    const real m = ami[i]*fei[i-1];
    fei[i] = 1.0/(aci[i] - m*api[i-1]);
    rrr[i] -= ami[i]*fei[i-1]*rrr[i-1];
    u [i] = -ami[i]*fei[i-1]* u [i-1];
  }
  rrr[n-1] *= fei[n-1];
  u [n-1] += api[n-1];
  u [n-1] *= fei[n-1];
  for (int i = n-2; i >= 0; i--)
  {
    rrr[i] = (rrr[i] - api[i]*rrr[i+1])*fei[i];
    u [i] = ( u [i] - api[i]* u [i+1])*fei[i];
  }
#endif

  /*      Form v · x/(1 + v · z). */

  const real fact = -1.0 *
    (rrr[0] + rrr[n-1]*ami[0]/gamma) /
    ( u [0] +  u [n-1]*ami[0]/gamma + 1.0);

  /*      Now get the solution vector x */

  for (int i = 0; i < n; i++)
    rrr[i] += fact*u[i];

  aci[ 0 ] = aci0;
  aci[n-1] = acin;
}
#endif
