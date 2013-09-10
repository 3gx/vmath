struct Force
{
  double *ax;
  double *ay;
  double *az;
  double *jx;
  double *jy;
  double *jz;
};

struct Predictor
{
  double *x;
  double *y;
  double *z;
  double *vx;
  double *vy;
  double *vz;
  double *mass;
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
    uniform       double gax[],
    uniform       double gay[],
    uniform       double gaz[],
    uniform       double gjx[],
    uniform       double gjy[],
    uniform       double gjz[])
{
  uniform const Predictor pred =  pred0[0];
  foreach (i = ibeg ... iend)
  {
    double ax, ay, az;
    double jx, jy, jz;
    ax = ay = az = 0.0;
    jx = jy = jz = 0.0;

    const double  xi = pred.x [i];
    const double  yi = pred.y [i];
    const double  zi = pred.z [i];
    const double vxi = pred.vx[i];
    const double vyi = pred.vy[i];
    const double vzi = pred.vz[i];

    for (uniform int j = jbeg; j < jend; j += dnj)
    {
      const double  xj = pred.x[j];
      const double  yj = pred.y[j];
      const double  zj = pred.z[j];
      const double vxj = pred.vx[j];
      const double vyj = pred.vy[j];
      const double vzj = pred.vz[j];
      const double  mj = pred.mass[j];

      const double dx = xj - xi;
      const double dy = yj - yi;
      const double dz = zj - zi;

      const double ds2 = dx*dx + dy*dy + dz*dz + eps2;

      const double inv_ds = rsqrt((float)ds2);
      const double  inv_ds2 = inv_ds  * inv_ds;
      const double minv_ds  = inv_ds  * mj;
      const double minv_ds3 = inv_ds2 * minv_ds;


      ax += minv_ds3 * dx;
      ay += minv_ds3 * dy;
      az += minv_ds3 * dz;

      const double dvx = vxj - vxi;
      const double dvy = vyj - vyi;
      const double dvz = vzj - vzi;
      const double rv  = dx*dvx + dy*dvy + dz*dvz;

      const double Jij = (-3.0) * (rv * inv_ds2 * minv_ds3);

      jx += minv_ds3*dvx + Jij*dx;
      jy += minv_ds3*dvy + Jij*dy;
      jz += minv_ds3*dvz + Jij*dz;
    }
    gax[i-ioffset] = ax;
    gay[i-ioffset] = ay;
    gaz[i-ioffset] = az;
    gjx[i-ioffset] = jx;
    gjy[i-ioffset] = jy;
    gjz[i-ioffset] = jz;
  }
}

export void calcforce_on_first_nact1(
    uniform const int nj,
    uniform const int i,
    uniform const double eps2,
    uniform const Predictor pred0[],
    uniform Force force0[])
{
  uniform const Predictor pred =  pred0[0];
  uniform       Force    force = force0[0];

  double ax, ay, az;
  double jx, jy, jz;
  ax = ay = az = 0.0;
  jx = jy = jz = 0.0;

  const double  xi = pred.x[i];
  const double  yi = pred.y[i];
  const double  zi = pred.z[i];
  const double vxi = pred.vx[i];
  const double vyi = pred.vy[i];
  const double vzi = pred.vz[i];

  foreach (j=0...nj)
  {
    const double  xj = pred. x[j];
    const double  yj = pred. y[j];
    const double  zj = pred. z[j];
    const double vxj = pred.vx[j];
    const double vyj = pred.vy[j];
    const double vzj = pred.vz[j];
    const double  mj = pred.mass[j];

    const double dx = xj - xi;
    const double dy = yj - yi;
    const double dz = zj - zi;

    const double ds2 = dx*dx + dy*dy + dz*dz + eps2;

    const double inv_ds = rsqrt((float)ds2);
    const double  inv_ds2 = inv_ds  * inv_ds;
    const double minv_ds  = inv_ds  * mj;
    const double minv_ds3 = inv_ds2 * minv_ds;


    ax += minv_ds3 * dx;
    ay += minv_ds3 * dy;
    az += minv_ds3 * dz;

    const double dvx = vxj - vxi;
    const double dvy = vyj - vyi;
    const double dvz = vzj - vzi;
    const double rv  = dx*dvx + dy*dvy + dz*dvz;

    const double Jij = (-3.0) * (rv * inv_ds2 * minv_ds3);

    jx += minv_ds3*dvx + Jij*dx;
    jy += minv_ds3*dvy + Jij*dy;
    jz += minv_ds3*dvz + Jij*dz;
  }

  force.ax[i] = reduce_add(ax);
  force.ay[i] = reduce_add(ay);
  force.az[i] = reduce_add(az);
  force.jx[i] = reduce_add(jx);
  force.jy[i] = reduce_add(jy);
  force.jz[i] = reduce_add(jz);
}

void body_body_force(
    double &ax, double &ay, double &az, 
    double &jx, double &jy, double &jz,
    const double xi, const double yi, const double zi,
    const double vxi, const double vyi, const double vzi,
    const double xj, const double yj, const double zj,
    const double vxj, const double vyj, const double vzj, const double mj,
    const double eps2)
{
  const double dx = xj - xi;
  const double dy = yj - yi;
  const double dz = zj - zi;

  const double ds2 = dx*dx + dy*dy + dz*dz + eps2;

  const double  inv_ds  = rsqrt(ds2);
  const double  inv_ds2 = inv_ds*inv_ds;
  const double minv_ds  = inv_ds  * mj;
  const double minv_ds3 = inv_ds2 * minv_ds;

  ax += minv_ds3 * dx;
  ay += minv_ds3 * dy;
  az += minv_ds3 * dz;
  
  const double dvx = vxj - vxi;
  const double dvy = vyj - vyi;
  const double dvz = vzj - vzi;
  const double rv  = dx*dvx + dy*dvy + dz*dvz;
  
  const double Jij = (-3.0) * (rv * inv_ds2 * minv_ds3);
  
  jx += minv_ds3*dvx + Jij*dx;
  jy += minv_ds3*dvy + Jij*dy;
  jz += minv_ds3*dvz + Jij*dz;
}

export void compute_forces(
    uniform const int     n,
    uniform const int     n1,
    uniform const int     n2,
    uniform const double eps2,
    uniform const double mass[],
    uniform const double posx[],
    uniform const double posy[],
    uniform const double posz[],
    uniform const double velx[],
    uniform const double vely[],
    uniform const double velz[],
    uniform       double accx[],
    uniform       double accy[],
    uniform       double accz[],
    uniform       double jrkx[],
    uniform       double jrky[],
    uniform       double jrkz[])
{
  foreach (i = n1 ... n2)
  {
    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;
    double jx = 0.0;
    double jy = 0.0;
    double jz = 0.0;
    const double  xi = posx[i];
    const double  yi = posy[i];
    const double  zi = posz[i];
    const double vxi = velx[i];
    const double vyi = vely[i];
    const double vzi = velz[i];

    for (uniform int j = 0; j < n; j++)
    {
      const double  xj = posx[j];
      const double  yj = posy[j];
      const double  zj = posz[j];
      const double vxj = velx[j];
      const double vyj = vely[j];
      const double vzj = velz[j];
      const double  mj = mass[j];
      body_body_force(ax,ay,az,jx,jy,jz,
          xi,yi,zi,vxi,vyi,vzi,
          xj,yj,zj,vxj,vyj,vzj,mj,eps2);
    }
    accx[i] = ax;
    accy[i] = ay;
    accz[i] = az;
    jrkx[i] = jx;
    jrky[i] = jy;
    jrkz[i] = jz;
  }
}

#undef double
