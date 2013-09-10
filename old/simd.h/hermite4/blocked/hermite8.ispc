struct Predictor
{
  double* pos[4];  /* x, y, z, mass */
  double* vel[3];
  double* acc[3];
  double* jrk[3];
};

struct Force
{
  double* acc[3];
  double* jrk[3];
  double* snp[3];
  double* crk[3];
};


export void calcforce_on_first_nact(
    uniform const int ibeg,
    uniform const int iend,
    uniform const int ioffset,
    uniform const int jbeg,
    uniform const int jend,
    uniform const int dnj,
    uniform const double eps2,
    uniform const Predictor pred0[],
    uniform       double gaccx[],
    uniform       double gaccy[],
    uniform       double gaccz[],
    uniform       double gjrkx[],
    uniform       double gjrky[],
    uniform       double gjrkz[],
    uniform       double gsnpx[],
    uniform       double gsnpy[],
    uniform       double gsnpz[],
    uniform       double gcrkx[],
    uniform       double gcrky[],
    uniform       double gcrkz[])
{
  uniform const Predictor pred =  pred0[0];
  foreach (i = ibeg ... iend)
  {
    double acc[3] = {0.0, 0.0, 0.0};
    double jrk[3] = {0.0, 0.0, 0.0};
    double snp[3] = {0.0, 0.0, 0.0};
    double crk[3] = {0.0, 0.0, 0.0};

    const double  xi = pred.pos[0][i];
    const double  yi = pred.pos[1][i];
    const double  zi = pred.pos[2][i];

    const double vxi = pred.vel[0][i];
    const double vyi = pred.vel[1][i];
    const double vzi = pred.vel[2][i];

    const double axi = pred.acc[0][i];
    const double ayi = pred.acc[1][i];
    const double azi = pred.acc[2][i];

    const double jxi = pred.jrk[0][i];
    const double jyi = pred.jrk[1][i];
    const double jzi = pred.jrk[2][i];

    for (uniform int j = jbeg; j < jend; j += dnj)
    {
      const double  xj = pred.pos[0][j];
      const double  yj = pred.pos[1][j];
      const double  zj = pred.pos[2][j];
      const double  mj = pred.pos[3][j];

      const double vxj = pred.vel[0][j];
      const double vyj = pred.vel[1][j];
      const double vzj = pred.vel[2][j];

      const double axj = pred.acc[0][j];
      const double ayj = pred.acc[1][j];
      const double azj = pred.acc[2][j];
    
      const double jxj = pred.jrk[0][j];
      const double jyj = pred.jrk[1][j];
      const double jzj = pred.jrk[2][j];

      const double dx = xj - xi;
      const double dy = yj - yi;
      const double dz = zj - zi;  

      const double dvx = vxj - vxi;
      const double dvy = vyj - vyi;
      const double dvz = vzj - vzi;

      const double dax = axj - axi;
      const double day = ayj - ayi;
      const double daz = azj - azi;
      
      const double djx = jxj - jxi;
      const double djy = jyj - jyi;
      const double djz = jzj - jzi;

      // 12

      const double r2   = dx*dx + dy*dy + dz*dz + eps2;

      const double drdv =  dx*dvx +  dy*dvy +  dz*dvz;
      const double drda =  dx*dax +  dy*day +  dz*daz;
      const double drdj =  dx*djx +  dy*djy +  dz*djz;
      const double dvdv = dvx*dvx + dvy*dvy + dvz*dvz;
      const double dvda = dvx*dax + dvy*day + dvz*daz;

      // 31

      double ri    = rsqrt(r2);
      ri = ((-0.5)*ri)*(r2*ri*ri + (-3.0));
#if 1
      ri = ((-0.5)*ri)*(r2*ri*ri + (-3.0));
#endif
      const double ri2   = ri*ri;
      const double mri3  = mj * ri2 * ri;
      const double alpha = drdv * ri2;
      const double beta  = (dvdv + drda) * ri2 + alpha * alpha;
      const double gamma = (3.0*dvda + drdj)*ri2 + alpha*(3.0*beta - 4.0*alpha*alpha);

      acc[0] += mri3 * dx;
      acc[1] += mri3 * dy;
      acc[2] += mri3 * dz;

      // 23

      const double c1 = (-3.0)*alpha;
      const double tmp1x = dvx + c1 * dx;
      const double tmp1y = dvy + c1 * dy;
      const double tmp1z = dvz + c1 * dz;

      jrk[0] += mri3 * tmp1x;
      jrk[1] += mri3 * tmp1y;
      jrk[2] += mri3 * tmp1z;

      //13

      const double c2 = (-6.0)*alpha;
      const double c3 = (-3.0)*beta;
      const double tmp2x = dax + c2 * tmp1x + c3 * dx;
      const double tmp2y = day + c2 * tmp1y + c3 * dy;
      const double tmp2z = daz + c2 * tmp1z + c3 * dz;
      snp[0] += mri3 * tmp2x;
      snp[1] += mri3 * tmp2y;
      snp[2] += mri3 * tmp2z;

      // 20
      
      const double c5 = (-9.0)*alpha;
      const double c6 = (-9.0)*beta;
      const double c7 = (-3.0)*gamma;
      const double tmp3x = djx + c5 * tmp2x + c6 * tmp1x + c7 * dx;
      const double tmp3y = djy + c5 * tmp2y + c6 * tmp1y + c7 * dy;
      const double tmp3z = djz + c5 * tmp2z + c6 * tmp1z + c7 * dz;

      crk[0] += mri3 * tmp3x;
      crk[1] += mri3 * tmp3y;
      crk[2] += mri3 * tmp3z;

      // 27
    } // 126
    gaccx[i-ioffset] = acc[0];
    gaccy[i-ioffset] = acc[1];
    gaccz[i-ioffset] = acc[2];
    gjrkx[i-ioffset] = jrk[0];
    gjrky[i-ioffset] = jrk[1];
    gjrkz[i-ioffset] = jrk[2];
    gsnpx[i-ioffset] = snp[0];
    gsnpy[i-ioffset] = snp[1];
    gsnpz[i-ioffset] = snp[2];
    gcrkx[i-ioffset] = crk[0];
    gcrky[i-ioffset] = crk[1];
    gcrkz[i-ioffset] = crk[2];
  }
}

