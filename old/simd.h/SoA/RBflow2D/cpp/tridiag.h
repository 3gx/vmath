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
