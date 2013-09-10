#include "defs_ispc.h"

uniform Data data;

export void
init(
    uniform real rc[], uniform real rm[], uniform real g2rc[], uniform real g2rm[],
    uniform real zz[], uniform real zm[], uniform real g3rc[], uniform real g3rm[],

    uniform real udx2rm[], uniform real udx2c[], uniform real udx2m[],
    uniform real usg2rc[], uniform real udx3c[], uniform real udx3m[],

    uniform int   jmv[], uniform int  jmc[], uniform int  jpv[], uniform int  jpc[], uniform int  jmhv[],
    uniform int   kmc[], uniform int  kpc[], uniform int  kmv[], uniform int  kpv[], uniform int  kup[], uniform int  kum[],

    uniform real ap2j[],  uniform real ac2j[],  uniform real am2j[],
    uniform real ap2je[], uniform real ac2je[], uniform real am2je[],
    uniform real ap3j[],  uniform real ac3j[],  uniform real am3j[],
    uniform real apscj[], uniform real acscj[], uniform real amscj[],

    uniform real ap3ck[],  uniform real ac3ck[],  uniform real am3ck[],
    uniform real ap3sk[],  uniform real ac3sk[],  uniform real am3sk[],
    uniform real ap3ssk[], uniform real ac3ssk[], uniform real am3ssk[],

    uniform real ak2[],
    uniform real amphj[], uniform real apphj[], uniform real acphj[],
    uniform real amphk[], uniform real acphk[], uniform real apphk[],

    uniform int  n2m,  uniform int n3m,
    uniform real ren, uniform real pec)
{
  data.rc   = rc;
  data.rm   = rm;
  data.g2rc = g2rc;
  data.g2rm = g2rm;

  data.zz   = zz;
  data.zm   = zm;
  data.g3rc = g3rc;
  data.g3rm = g3rm;


  data.udx2rm = udx2rm;
  data.udx2c  = udx2c;
  data.udx2m  = udx2m;

  data.usg2rc = usg2rc;
  data.udx3c  = udx3c;
  data.udx3m  = udx3m;


  data.jmv  = jmv;
  data.jmc  = jmc;
  data.jpv  = jpv;
  data.jmhv = jmhv;
  data.kmc  = kmc;
  data.kpc  = kpc;
  data.kmv  = kmv;
  data.kpv  = kpv;  
  data.kup  = kup;
  data.kum  = kum;


  data.ap2j  = ap2j;
  data.ac2j  = ac2j;
  data.am2j  = am2j;
  data.ap2je = ap2je;
  data.ac2je = ac2je;
  data.am2je = am2je;
  data.ap3j  = ap3j;
  data.ac3j  = ac3j;
  data.am3j  = am3j;
  data.apscj = apscj;
  data.acscj = acscj;
  data.amscj = amscj;


  data.ap3ck  = ap3ck;
  data.ac3ck  = ac3ck;
  data.am3ck  = am3ck;
  data.ap3sk  = ap3sk;
  data.ac3sk  = ac3sk;
  data.am3sk  = am3sk;
  data.ap3ssk = ap3ssk;
  data.ac3ssk = ac3ssk;
  data.am3ssk = am3ssk;

  data.ak2   = ak2;
  data.amphj = amphj;
  data.acphj = acphj;
  data.apphj = apphj;
  data.amphk = amphk;
  data.acphk = acphk;
  data.apphk = apphk;

  data.n2m = n2m;
  data.n3m = n3m;
  data.ren = ren;
  data.pec = pec;
}

export void 
set_algaro(uniform real al, uniform real ga, uniform real ro)
{
  data.al = al;
  data.ga = ga;
  data.ro = ro;
}

export void 
set_data( 
    uniform real dt, uniform real beta, 
    uniform real denbs[], uniform real denbn[],
    uniform real q2[],  uniform real q3[],
    uniform real pr[],  uniform real dens[], uniform real hro[], 
    uniform real rhs1[], uniform real rhs2[], uniform real rhs3[],
    uniform real ru1[], uniform real ru2[],  uniform real ru3[], uniform real ruro[],
    uniform real qcap[], uniform real dph[])
{
  data.dt   = dt;
  data.beta = beta;

  data.denbs = denbs;
  data.denbn = denbn;

  data.q2   = (real (*)[M2])q2;
  data.q3   = (real (*)[M2])q3;
  data.pr   = (real (*)[M2])pr;
  data.dens = (real (*)[M2])dens;
  data.hro  = (real (*)[M2])hro;
  data.rhs1 = (real (*)[M2])rhs1;
  data.rhs2 = (real (*)[M2])rhs2;
  data.rhs3 = (real (*)[M2])rhs3;
  data.ru1  = (real (*)[M2])ru1;
  data.ru2  = (real (*)[M2])ru2;
  data.ru3  = (real (*)[M2])ru3;
  data.ruro = (real (*)[M2])ruro;
  data.qcap = (real (*)[M2])qcap;
  data.dph  = (real (*)[M2+16])dph;
}
