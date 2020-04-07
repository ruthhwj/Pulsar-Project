//START REGION RELEASE
#include <ctype.h>

//START REGION DEVELOP
//START REGION RELEASE
extern void sla_altaz(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

void csla_altaz(double ha, double dec, double phi, double *az, double *azd, double *azdd, double *el, double *eld, double *eldd, double *pa, double *pad, double *padd)
{
  sla_altaz(&ha, &dec, &phi, az, azd, azdd, el, eld, eldd, pa, pad, padd);
}


//START REGION DEVELOP


extern void sla_cldj(int *iy, int *im, int *id, double *djm, int *j);

void csla_cldj(int iy, int im, int id, double *djm, int *j)
{
  sla_cldj(&iy, &im, &id, djm, j);
}

//START REGION DEVELOP
//START REGION RELEASE

extern void sla_dcc2s(double *v, double *ra, double *dec);

void csla_dcc2s(double v[3], double *a, double *b)
{
  sla_dcc2s(v, a, b);
}

//START REGION DEVELOP
//START REGION RELEASE

extern void sla_dcs2c(double *ra, double *dec, double *v);

void csla_dcs2c(double a, double b, double v[3])
{
  sla_dcs2c(&a, &b, v);
}

//START REGION DEVELOP

extern void sla_de2h(double *ha, double *dec, double *phi, double *az, double *el);

void csla_de2h(double ha, double dec, double phi, double *az, double *el)
{
  sla_de2h(&ha, &dec, &phi, az, el);
}


//START REGION DEVELOP
//START REGION RELEASE

extern void sla_dmxv(double *dm, double *va, double *vb);

void csla_dmxv(double dm[3][3], double va[3], double vb[3])
{
  sla_dmxv(&(dm[0][0]), va, vb);
}

//START REGION DEVELOP
//START REGION RELEASE

extern double sla_dranrm(double *);

double csla_dranrm(double angle)
{
  return sla_dranrm(&angle);
}

//START REGION DEVELOP

extern double sla_dsep(double *, double *, double *, double*);

double csla_dsep(double a1, double b1, double a2, double b2)
{
  return sla_dsep(&a1, &b1, &a2, &b2);
}



extern double sla_dtt(double *mjd);

double csla_dtt(double dju)
{
  return sla_dtt(&dju);
}


//START REGION DEVELOP
//START REGION RELEASE

extern double sla_dvdv(double *va, double *vb);

double csla_dvdv(double va[3], double vb[3])
{
  return sla_dvdv(va, vb);
}

//START REGION DEVELOP
//START REGION RELEASE

extern double sla_dvn(double *v, double *uv, double *vm);

double csla_dvn(double v[3], double uv[3], double *vm)
{
  return sla_dvn(v, uv, vm);
}

//START REGION DEVELOP

extern void sla_ecleq(double *, double *, double *, double *, double *);

void csla_ecleq(double dl, double db, double mjd, double *dr, double *dd)
{
  sla_ecleq(&dl, &db, &mjd, dr, dd);
}

//START REGION DEVELOP
//START REGION RELEASE

extern double sla_epb(double *mjd);

double csla_epb(double date)
{
  return sla_epb(&date);
}

//START REGION RELEASE

extern double sla_epj(double *mjd);

double csla_epj(double date)
{
  return sla_epj(&date);
}

//START REGION DEVELOP

extern void sla_eqecl(double *, double *, double *, double *, double *);

void csla_eqecl(double dr, double dd, double mjd, double *dl, double *db)
{
  sla_eqecl(&dr, &dd, &mjd, dl, db);
}

extern void sla_eqgal(double *dr, double *dd, double *dl, double *db);

void csla_eqgal(double dr, double dd, double *dl, double *db)
{
  sla_eqgal(&dr, &dd, dl, db);
}

//START REGION DEVELOP
//START REGION RELEASE

extern void sla_evp(double *tdb, double *ep, double *dvb, double *dpb, double *dvh, double *dph);

void csla_evp(double date, double deqx, double dvb[3], double dpb[3], double dvh[3], double dph[3])
{
  sla_evp(&date, &deqx, dvb, dpb, dvh, dph);
}

//START REGION DEVELOP


extern void sla_galeq(double *dl, double *db, double *dr, double *dd);

void csla_galeq(double dl, double db, double *dr, double *dd)
{
  sla_galeq(&dl, &db, dr, dd);
}

extern void sla_geoc(double *p, double *h, double *r, double *z);

void csla_geoc(double p, double h, double *r, double *z)
{
  sla_geoc(&p, &h, r, z);
}

//START REGION DEVELOP
//START REGION RELEASE

extern double sla_gmst(double *mjd);

double csla_gmst(double ut1)
{
  return sla_gmst(&ut1);
}

//START REGION DEVELOP

extern void sla_map(double *rm, double *dm, double *pr, double *pd, double *px, double *rv, double *eq, double *date, double *ra, double *da);

void csla_map(double rm, double dm, double pr, double pd, double px, double rv, double eq, double date, double *ra, double *da)
{
  sla_map(&rm, &dm, &pr, &pd, &px, &rv, &eq, &date, ra, da);
}

//START REGION DEVELOP
//START REGION RELEASE

extern void sla_nut(double *date, double *rmatn);

void csla_nut(double date, double rmatn[3][3])
{
  sla_nut(&date, &(rmatn[0][0]));
}

//START REGION DEVELOP
//START REGION RELEASE

extern double sla_pa(double *HA, double *DEC, double *PHI);

double csla_pa(double ha, double dec, double phi)
{
  return sla_pa(&ha, &dec, &phi);
}

//START REGION DEVELOP
//START REGION RELEASE

extern void sla_prebn(double *ep0, double *ep1, double *rmatp);

void csla_prebn(double bep0, double bep1, double rmatp[3][3])
{
  sla_prebn(&bep0, &bep1, &(rmatp[0][0]));
}

//START REGION DEVELOP
//START REGION RELEASE

extern void sla_prec(double *ep0, double *ep1, double *rmatp);

void csla_prec(double ep0, double ep1, double rmatp[3][3])
{
  sla_prec(&ep0, &ep1, &(rmatp[0][0]));
}

//START REGION DEVELOP

extern void sla_preces(char *system, double *ep0, double *ep1, double *ra, double *dec);

void csla_preces(char *system, double ep0, double ep1, double *ra, double *dec)
{
  /* Cannot get it to work with string, so translate the function instead */
  //    sla_preces(system, &ep0, &ep1, ra, dec);

  double pm[3][3], v1[3], v2[3];

  /* Validate sys */
  if(( toupper( system[0] ) != 'F' ) ||( toupper( system[1] ) != 'K' ) ||( system[2] != '4' && system[2] != '5' ) ) {
    *ra = -99.0;          /* Error */
    *dec = -99.0;
  } else {
    /* Generate appropriate precession matrix */
    if( system[2] == '4' )
      csla_prebn( ep0, ep1, pm );
    else
      csla_prec( ep0, ep1, pm );

   /* Convert RA,Dec to x,y,z */
    csla_dcs2c( *ra, *dec, v1 );

   /* Precess */
    csla_dmxv( pm, v1, v2 );

   /* Back to RA,Dec */
    csla_dcc2s( v2, ra, dec );
    *ra = csla_dranrm( *ra );
  }
}

extern void sla_precl(double *ep0, double *ep1, double *rmatp);

void csla_precl(double ep0, double ep1, double rmatp[3][3])
{
  sla_precl(&ep0, &ep1, &(rmatp[0][0]));
}

//START REGION DEVELOP
//START REGION RELEASE

extern void sla_prenut(double *epoch, double *date, double *rmatpn);

void csla_prenut(double ep0, double date, double rmatpn[3][3])
{
  sla_prenut(&ep0, &date, &(rmatpn[0][0]));
}

//START REGION DEVELOP

extern void sla_subet(double *rc, double *dc, double *eq, double *rm, double *dm);

void csla_subet(double rc, double dc, double eq, double *rm, double *dm)
{
  sla_subet(&rc, &dc, &eq, rm, dm);
}





//START REGION DEVELOP
