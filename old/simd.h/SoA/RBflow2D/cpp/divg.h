#pragma once

template<typename REAL>
inline void divg_ij(
    const REAL usdtal,
    const int jc, const int kc,
    const int jp, const int kp)
{
  const REAL dqcap = 
    (LDU(q2(jp,kc)) - LDA(q2(jc,kc))) * LDA  (udx2m(jc)) + 
    (LDA(q3(jc,kp)) - LDA(q3(jc,kc))) * REAL(udx3m(kc));

  REF(dph(jc,kc)) = dqcap*usdtal;
}

void divg()
{
  assert(n2m % vreal::WIDTH == 0);
  
  const real usdtal = 1.0/(dt*al);

#pragma omp parallel for
  for (int kc = 0; kc < n3m; kc++)
    foreach_warp(jc, 0, n2m)
      if (jc < vreal::WIDTH || jc >= n2m - vreal::WIDTH)
        foreach_lane(ch)
        {
          const int jp = jpv(jc+ch) - 1;
          divg_ij<sreal>(usdtal, jc+ch, kc, jp, kc+1);
        }
      else
        divg_ij<vreal>(usdtal, jc, kc, jc+1, kc+1);
}


