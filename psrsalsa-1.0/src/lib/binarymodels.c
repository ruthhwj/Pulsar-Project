// The following functions are adaptations of those implemented in tempo2

//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
*    This file is part of TEMPO2. 
* 
*    TEMPO2 is free software: you can redistribute it and/or modify 
*    it under the terms of the GNU General Public License as published by 
*    the Free Software Foundation, either version 3 of the License, or 
*    (at your option) any later version. 
*    TEMPO2 is distributed in the hope that it will be useful, 
*    but WITHOUT ANY WARRANTY; without even the implied warranty of 
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
*    GNU General Public License for more details. 
*    You should have received a copy of the GNU General Public License 
*    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
*/

/*
*    If you use TEMPO2 then please acknowledge it by citing 
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model.
*/

#include <math.h>

#include <psrsalsa.h>

#define SecPerDay           86400.0         /* 24.0*3600.0 */

/* Based on bnrybt.f, only considering one orbit */
// bbat = mjd as a bat (in days)
// t0 = Epoch of periastron (MJD)
// pb = Orbital period (days)
// ecc = Orbital eccentricity
// pbdot = The first time derivative of binary period, same unit as tempo2 (s/s?)
// a1 = Projected semimajor axis of orbit (lt-s)
// om = Longitude of periastron (deg)
// omdot = Periastron advance (deg/yr) 
// a1dot = Rate of change of semimajor axis (lt-s / s) 
// gamma = Post-Keplerian "gamma" term (s) 
//
// It returns the time delay caused by the orbit, i.e. the expected residual.

double BTmodel(long double bbat, binary_def companion)
{
  double torb;
  double tt0;
  double orbits;
  double xpbdot;
  double edot;
  double asini;
  double xdot;
  double omega;
  int    norbits;
  double phase;
  double ep,dep,bige,tt,som,com;
  double alpha,beta,sbe,cbe,q,r,s;

  long double t0 = companion.t0;
  long double pb = companion.pb;
  long double ecc = companion.ecc;
  long double a1dot = companion.a1dot;
  long double a1 = companion.a1;
  long double om = companion.om;
  long double omdot = companion.omdot;
  long double pbdot = companion.pbdot;
  long double gamma = companion.gamma;

  tt0 = (bbat - t0)*86400.0;

  pb     *= 86400.0;
  edot   = 0.0;
  ecc    += edot*tt0;

  if (ecc < 0.0 || ecc > 1.0) {
    printerror(0, "ERROR tempo3: problem with eccentricity = %Lg\n", ecc);
    exit(0);
  }

  xpbdot = 0.0;
  xdot = a1dot;
  asini  = a1 + xdot*tt0;
  om *= M_PI/180.0;
  omega  = om + omdot*tt0/(86400.0 * 365.25)/(180.0/M_PI);

  torb = 0.0;
  /* Should ct be the barycentric arrival time? -- check bnrybt.f */
  orbits = tt0/pb - 0.5*(pbdot+xpbdot)*pow(tt0/pb,2); 
  norbits = (int)orbits;
  if (orbits < 0.0) norbits--;
  
  phase = 2.0*M_PI * (orbits-norbits);

  /* Using Pat Wallace's method of solving Kepler's equation -- code based on bnrybt.f */
  ep = phase + ecc*sin(phase)*(1.0+ecc*cos(phase));

  /* This line is wrong in the original tempo: should be inside the do loop */
  /*  denom = 1.0 - ecc*cos(ep);*/
  
  do {
    dep = (phase - (ep-ecc*sin(ep)))/(1.0 - ecc*cos(ep));
    ep += dep;
  } while (fabs(dep) > 1.0e-12);
  bige = ep;

  tt = 1.0-ecc*ecc;
  som = sin(omega);
  com = cos(omega);

  alpha = asini*som;
  beta = asini*com*sqrt(tt);
  sbe = sin(bige);
  cbe = cos(bige);
  q = alpha * (cbe-ecc) + (beta+gamma)*sbe;
  r = -alpha*sbe + beta*cbe;
  s = 1.0/(1.0-ecc*cbe);

  torb = -q+(2*M_PI/pb)*q*r*s + torb;

  return torb;

  /*
  if (param==-1) return torb;

  if (param==param_pb)
    return -2.0*M_PI*r*s/pb*86400.0*tt0/(86400.0*pb) * 86400.0;  // fctn(12+j) 
  else if (param==param_a1)
    return (som*(cbe-ecc) + com*sbe*sqrt(tt));                // fctn(9+j) 
  else if (param==param_ecc)
    return -(alpha*(1.0+sbe*sbe-ecc*cbe)*tt - beta*(cbe-ecc)*sbe)*s/tt; // fctn(10+j) 
  else if (param==param_om)
    return asini*(com*(cbe-ecc) - som*sqrt(tt)*sbe);          // fctn(13+j) 
  else if (param==param_t0)
    return -2.0*M_PI/pb*r*s*86400.0;                           // fctn(11+j) 
  else if (param==param_pbdot)
    return 0.5*(-2.0*M_PI*r*s/pb*86400.0*tt0/(86400.0*pb))*tt0; // fctn(18+j) 
  else if (param==param_a1dot)
    return (som*(cbe-ecc) + com*sbe*sqrt(tt))*tt0;            // fctn(24+j) 
  else if (param==param_omdot)
    return asini*(com*(cbe-ecc) - som*sqrt(tt)*sbe)*tt0;      // fctn(14+j) 
  else if (param==param_edot)                            
    return (-(alpha*(1.0+sbe*sbe-ecc*cbe)*tt - beta*(cbe-ecc)*sbe)*s/tt)*tt0; // fctn(25+j) 
  else if (param==param_gamma) 
    return sbe;                                               // fctn(15+j) 
  return 0.0;
*/
}


// This is a direct implementation of tempo2's T2model.C


// bbat = mjd as a bat (in days)
// t0 = Epoch of periastron (MJD)
// pb = Orbital period (days)
// ecc = Orbital eccentricity
// pbdot = The first time derivative of binary period, same unit as tempo2 (s/s?)
// a1 = Projected semimajor axis of orbit (lt-s)
// om = Longitude of periastron (deg)
// omdot = Periastron advance (deg/yr) 
// a1dot = Rate of change of semimajor axis (lt-s / s) 
// gamma = Post-Keplerian "gamma" term (s) 
// eps1= ELL1 binary model parameter 1 
// eps2= ELL1 binary model parameter 2 
// tasc= Time of ascending node
// shapmax = ....
// kom = .... in deg
// kin = .... in deg
// DIDN'T IMPLEMENT: bpjep related parameters, code is commented out. These are jumps in binary parameters. 
//                   reading in of many other (post-Keplerian) parameters, but used in calculations
//                   beta-prime parameters (also not implemented in tempo2 for this model)
//                   mtot calculation not implemented
//                   Use of the KopeikinTerms not implemented
//                   Shapiro delay calculation not implemented: h3, stig, h4

// nCompanion = I assume it is the number of companions
// Not sure what the variable arr is, but it might be to select an individual companion only. I think a parameter number we don't use.
// It returns the time delay caused by the orbit, i.e. the expected residual.

double T2model(long double bbat, int nCompanion, binary_def *companion)
{
  double an;
  double pb,omdot;
  double rad2deg = 180.0/M_PI;
  double m2,tt0,t0,x,ecc,er,xdot,edot,dr,dth,eth,ct;
  double pbdot,xpbdot,phase,u,gamma;
  double orbits;
  int    norbits;
  double cu,onemecu=0,cae,sae,ae,omega,omz,sw,cw,alpha,beta,bg,dre,drep,drepp,
    anhat,su=0;
  double sqr1me2,cume,brace,si,dlogbr,ds,da,a0,b0,d2bar,torb;
  double eps1,eps2,eps1dot,eps2dot,t0asc;
  double shapmax,sdds;
  int    com,com1,com2;
  int    allTerms=1;            /* = 0 to emulate BT model */
  double dpara;
  double pmra,pmdec;
  //  double sin_omega,cos_omega
  double ki;
  double mtot,xomdot,afac,kom; //,m1

  double daop;// DAOP is the time delay due to annual orbital
              // parallax. daop is the aop distance.
  long double DAOP=0.0L, DSR=0.0L;
  long double DOP=0.0L; // Orbital parallax time delay

  /*
  double SUNMASS = 4.925490947e-6;
  double ceps1,ceps2;
    double xk,cshapmax;
  double csigma,ce,cx,comega,cgamma,cdth,cm2,csi, ckom, ckin;

  long double DK011,DK012, DK021, DK022, DK031,DK032, DK041,DK042,C,S;
  long double DK013, DK014, DK023, DK024, DK033, DK034, DK043, DK044;

  long double h3,h4,stig;
  long double lgf,TrueAnom;
  long double lsc, fs;
  int nharm=4;
  int mode = -1; // See ELL1Hmodel.C
  */

  torb = 0.0;

  //  if (param==-1) 
  //    {
  // I assume this means: include all companions
  com1 = 0;
  com2 = nCompanion;
      //    }
      //  else
      //    {
      //      com1 = arr;
      //      com2 = arr+1;
      //    }

  //    printf("Number of companions = %d %d\n",com1,com2);
  
  for (com = com1; com < com2; com++)
    {      
      /* Obtain Keplerian parameters */
      //      getKeplerian(&psr[p],com,&pb,&t0,&ecc,&omz,&x,&eps1,&eps2,&t0asc,
      //                   &shapmax,&kom,&ki);
      pb  = companion[com].pb*SecPerDay;
      t0  = companion[com].t0;
      ecc = companion[com].ecc; 
      omz = companion[com].om;
      x   = companion[com].a1;
      eps1 = companion[com].eps1;
      eps2 = companion[com].eps2;
      t0asc    = companion[com].tasc;
      //      shapmax = companion[com].shapmax;
      //      kom = companion[com].kom*M_PI/180.0;
      //      ki      = companion[com].kin*M_PI/180.0;
      shapmax = kom = ki = 0;

      /* Now add in the jumps */
      //      addKeplerianJumps(&psr[p],ipos,&torb,&x,&ecc,&omz,&pb);
      /*
       * Following the BTJ model, this function adds jumps to the Keplerian parameters
       * at a specified epoch
       *
       */
      /*
      void addKeplerianJumps(pulsar *psr,int ipos,double *torb,double *x,double *ecc,
			     double *omz,double *pb)
      {
	int i;
	
	for (i=0;i<psr->param[param_bpjep].aSize;i++)
	  {
	    if (psr->param[param_bpjep].paramSet[i]==1 && 
		psr->obsn[ipos].bbat > psr->param[param_bpjep].val[i])
	      {
		torb = torb - (double)(psr->param[param_bpjph].val[i]
					 / psr->param[param_f].val[0]);  
		x    = x    + (double)psr->param[param_bpja1].val[i];
		ecc  = ecc  + (double)psr->param[param_bpjec].val[i];
		omz  = omz  + (double)psr->param[param_bpjom].val[i];    
		pb   = pb   + (double)psr->param[param_bpjpb].val[i]*SECDAY; 
	      }
	  } 
      }
      */

      /* Parameters derived from the Keplerian parameters */
      //      deriveKeplerian(pb,kom,&an,&sin_omega,&cos_omega);
      //      void deriveKeplerian(double pb,double kom,double *an,double *sin_omega,
      //                     double *cos_omega){
      an    = 2.0*M_PI/pb;
      //      sin_omega = sin(kom);
      //      cos_omega = cos(kom);


      /* Obtain post-Keplerian parameters */
      //      getPostKeplerian(&psr[p],com,an,&si,&m2,&mtot,&omdot,&gamma,&xdot,&xpbdot,
      //                       &pbdot,&edot,&pmra,&pmdec,&dpara,&dr,&dth,&a0,&b0,
      //                       &xomdot,&afac,&eps1dot,&eps2dot,&daop);
      // IGNORED MOST FOR NOW
      si=m2=mtot=xpbdot=edot=pmra=pmdec=dpara=dr=dth=a0=b0=xomdot=afac=eps1dot=eps2dot=daop = 0;

      omdot = companion[com].omdot*M_PI/(180.0*365.25*SecPerDay*an); 
      pbdot = companion[com].pbdot; 
      xdot  = companion[com].a1dot; 
      gamma = companion[com].gamma; 

      /*
  double SUNMASS = 4.925490947e-6;
  double rad2deg = 180.0/M_PI;
  //double pxConv = 1.74532925199432958E-2/3600.0e3;//converts mas to rad
  double pxConv = M_PI/180.0/3600*1e-3; // converts mas to rad
  double daopConv = 3.08568025e16;//pc in m

  //  logdbg("Going to get parameters");
  *si      = getParameter(psr,param_sini,com);
  if (*si > 1.0)
    {
      displayMsg(1,"BIN1","SIN I > 1.0, setting to 1: should probably use DDS model","",psr[0].noWarnings);
      *si = 1.0;
      psr[0].param[param_sini].val[0] = 1.0;
    }
  if (*si < -1.0)
    {
      displayMsg(1,"BIN1","SIN I < -1.0, setting to -1: should probably use DDS model","",psr[0].noWarnings);
      *si = -1.0;
      psr[0].param[param_sini].val[0] = -1.0;
    }
  *m2      = getParameter(psr,param_m2,com)*SUNMASS;
  *mtot    = getParameter(psr,param_mtot,com)*SUNMASS;
  *xpbdot  = getParameter(psr,param_xpbdot,com);
  *edot    = getParameter(psr,param_edot,com);
  *pmra    = getParameter(psr,param_pmra,com)
    * M_PI/(180.0*3600.0e3)/(365.25*86400.0);
  *pmdec   = getParameter(psr,param_pmdec,com)
    * M_PI/(180.0*3600.0e3)/(365.25*86400.0);
  *dpara   = getParameter(psr,param_px,com)*pxConv;
  *dr      = getParameter(psr,param_dr,com);
  *dth     = getParameter(psr,param_dth,com);
  *a0      = getParameter(psr,param_a0,com);
  *b0      = getParameter(psr,param_b0,com);
  *xomdot  = getParameter(psr,param_xomdot,com)/(an*rad2deg*365.25*86400.0);
  *afac    = getParameter(psr,param_afac,com);
  *eps1dot = getParameter(psr,param_eps1dot,com);
  *eps2dot = getParameter(psr,param_eps2dot,com);
  *daop    = getParameter(psr,param_daop,com)*1e-3/pxConv;
       */
      




      /* If the beta-prime parameters are set then omdot, gamma, si, dr, er, 
         dth, a0 and b0 can be calculated from the beta-prime values
       * - see papers of Taylor, Wolszczan, Damour & Weisberg 1992 - Nature
       *                 Damour & Taylor 1992 - Phys Rev D
       */
      /*      if (psr[p].param[param_bp].paramSet[0]==1 && 
          psr[p].param[param_bpp].paramSet[0]==1)
        {
          //	  useBeta(psr[p]); 
          printf("Beta model not implemented yet\n");
        }
      */

      /* If general relativity is assummed to be correct and 
       * the total system mass has been determined then the values 
       * of dr,dth,er, eth, sini, gamma and pbdot can be calculated
       */
      if (mtot != 0)
        {
	  printerror(0, "ERROR T2model: Use of the total mass mtot parameter is not supported in the code (yet?).");
	  exit(0);
	  //          calcGR(mtot,m2,x,ecc,an,afac,(double)psr[p].param[param_f].val[0],
	  //                 &dr,&dth,&er,&eth,&xk,&si,&gamma,&pbdot,&a0,&b0);
	  //          omdot   = xomdot + xk; 
        }

      /* Derive parameters from the post-Keplerian parameters */
      //      derivePostKeplerian(mtot,m2,dr,dth,ecc,&m1,&er,&eth);
      //      void derivePostKeplerian(double mtot,double m2,double dr,double dth,
      //                         double ecc,double *m1,double *er,double *eth)
      //      m1  = mtot - m2;  /* Pulsar mass */
      er  = ecc*(1.0+dr);
      eth = ecc*(1.0+dth);




      /* Obtain delta T */
      ct  = bbat;      
      if (t0 != 0)        
        tt0 = (ct-t0)*SecPerDay;
      else if( t0asc != 0) 
        tt0 = (ct-t0asc)*SecPerDay;
      else {
        printf("ERROR T2model: T0 or TASC needs to be set in the parameter file\n");
        exit(1);
      }
      /* Update parameters with their time derivatives */
      //      updateParameters(edot,xdot,eps1dot,eps2dot,tt0,&ecc,&x,&eps1,&eps2);
      //      void updateParameters(double edot,double xdot,double eps1dot,double eps2dot,
      //                      double tt0,double *ecc,double *x,double *eps1,
      //                      double *eps2){
      (ecc)  += edot*tt0;
      (x)        += xdot*tt0;
      (eps1) += eps1dot*tt0;
      (eps2) += eps2dot*tt0;



      /* Do some checks */
      if (ecc < 0.0 || ecc > 1.0)
        {
	  printerror(0, "ERROR T2model: Eccentricity of orbit < 0 or > 1, quiting.");
	  exit(0);
	  //          psr[p].param[param_ecc].val[com]=0.0;
	  //          ecc = 0.0;
        }

      /* Obtain number of orbits in tt0 */
      orbits  = tt0/pb - 0.5*(pbdot+xpbdot)*pow(tt0/pb,2);
      norbits = (int)orbits;
      if (orbits<0.0) norbits--;
      
      /* Obtain phase of orbit */
      phase=2.0*M_PI*(orbits-norbits);
      
      if (companion[com].ecc_defined)
        {
          /* Compute eccentric anomaly u by iterating Kepler's equation */	 
	  u = getEccentricAnomaly(phase, ecc);

	  /*  DD equations 17b and 17c */
	  su=sin(u); cu=cos(u);
	  onemecu=1.0-ecc*cu;
	  cae=(cu-ecc)/onemecu;                         /* Equation 17b */
	  sae=sqrt(1.0-pow(ecc,2))*su/onemecu;          /* Equation 17c */
	  ae=atan2(sae,cae);  // True anomaly
	  if(ae<0.0) ae=ae+2.0*M_PI;
	  ae=2.0*M_PI*orbits + ae - phase;
	  omega=omz/rad2deg + omdot*ae;
	  sw=sin(omega);
	  cw=cos(omega);
	  /* DD equations 26, 27, 57: */
	  sqr1me2=sqrt(1.0-pow(ecc,2));
	  cume=cu-ecc;
	  /* Update parameters due to proper motion - Kopeikin 1996 */
	  /* And annual-orbital and orbital parallax - Kopeikin 1995 */
	  if (ki!=0 && kom !=0 && (pmra != 0 || pmdec != 0))
	    {
	      printerror(0, "ERROR T2model: Use of the KopeikinTerms in combination of proper motion is not supported in the code (yet?).");
	      exit(0);
	      /*
		KopeikinTerms(&psr[p],ipos,ki,pmra,sin_omega,pmdec,cos_omega,tt0,
		dpara,daop,si,&x,&DK011,&DK012,&DK021,&DK022,
		&DK031,&DK032,&DK041,&DK042,&DK013,&DK014,&DK023,
		&DK024,&DK033,&DK034,&DK043,&DK044);
		//              logdbg("did KopeikinTerms");
		C = (longdouble)(cw*(cu-er)-sqrt(1.0-pow(eth,2.0))*sw*su);
		S = (longdouble)(sw*(cu-er)+cw*sqrt(1.0-pow(eth,2.0))*su);
		DAOP = (DK011+DK012)*C-(DK021+DK022)*S;
		DSR = (DK031+DK032)*C+(DK041+DK042)*S;
		DOP = dpara/AULTSC/2.0*pow(x,2.0)*
		( pow( sin( ki ), -2.0 ) - 0.5 + 0.5 * pow( ecc, 2.0 ) *
		( 1 + pow( sw, 2.0 ) - 3.0 / pow( sin( ki ), 2.0 ) ) -
		2.0 * ecc * ( pow( sin( ki ), -2.0 ) - pow( sw, 2.0 ) ) * 
		cume - sqr1me2 * 2 * sw * cw * su * cume + 0.5 * 
		( cos( 2.0 * omega ) + pow( ecc, 2.0 ) * 
		( pow( sin( ki ), -2.0 ) + pow( cu, 2.0 ) ) ) *
		cos( 2.0 * u ) );
		
		//	      printf("T2model: %g DAOP = %g DSR = %g DOP = %g\n",(double)psr[p].obsn[ipos].bbat,(double)DAOP,(double)DSR,(double)DOP);
		
		//                logdbg("DAOP is %g and DSR is %g\n", (double)DAOP, (double)DSR);
		//                logdbg("DAOP is %g, DK011 and DK021 are %f and %f\n",
		//                       (double)DAOP,(double)DK011,(double)DK021);
		*/
	    }
      
      
	  if (shapmax != 0)  /* DDS model */
	    {
	      sdds  = 1.0 - exp(-1.0*shapmax);
	      brace = onemecu-sdds*(sw*cume+sqr1me2*cw*su);
	    }
	  else
	    brace=onemecu-si*(sw*cume+sqr1me2*cw*su);  
	  
	  da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + 
                                               ecc*cw); /* Equation 27 */
	  
	  /* DD equations 46 to 51 */	  
	  alpha=x*sw;                                   /* Equation 46  */
	  beta=x*sqrt(1-pow(eth,2))*cw;                 /* Equation 47  */
	  bg=beta+gamma;
	  dre=alpha*(cu-er) + bg*su;                    /* Equation 48  */
	  drep=-alpha*su + bg*cu;                       /* Equation 49  */
	  drepp=-alpha*cu - bg*su;                      /* Equation 50  */
	  anhat=an/onemecu;                             /* Equation 51  */
	  
	  dlogbr=log(brace);
	  ds=-2*m2*dlogbr;        /* Equation 26 */
	  
	}
      else if (companion[com].eps1_defined)  /* ELL1 model */
        {
          dre  = x*(sin(phase)-0.5*(eps1*cos(2.0*phase)-eps2*sin(2.0*phase)));
          drep = x*cos(phase);
          drepp=-x*sin(phase);
	  //          logdbg("going to Kopeikin");
          /* Update parameters due to proper motion - Kopeikin 1996 */
	  if (ki!=0 && kom !=0 && (pmra != 0 || pmdec != 0))
	    {
	      printerror(0, "Use of the KopeikinTerms in combination of proper motion is not supported in the code (yet?).");
	      exit(0);
	      /*
              S = (sin(phase)-0.5*(eps1*cos(2.0*phase)-eps2*sin(2.0*phase)));
              C = cos(phase)+0.5*(eps2*cos(2.0*phase)+eps1*sin(2.0*phase));
              KopeikinTerms(&psr[p],ipos,ki,pmra,sin_omega,pmdec,cos_omega,tt0,
                            dpara,daop,si,&x,&DK011,&DK012,&DK021,&DK022,
                            &DK031,&DK032,&DK041,&DK042,&DK013,&DK014,&DK023,
                            &DK024,&DK033,&DK034,&DK043,&DK044);
              DAOP = (DK011+DK012)*C-(DK021+DK022)*S;
              DSR = (DK031+DK032)*C+(DK041+DK042)*S;
	      */
            }
	  
          brace=1-si*sin(phase);
          da=a0*sin(phase)+b0*cos(phase);  
	  
          anhat = an; ecc = 0.0;
	  
	  
	  /* Shapiro delay */
	  /*
	  if ( psr[p].param[param_h3].paramSet[0] * psr[p].param[param_stig].paramSet[0] == 1 
	       || psr[p].param[param_h3].paramSet[0] * psr[p].param[param_h4].paramSet[0] == 1){
	    
	    //	    printf("Using the Friere & Wex formalism for the Shapiro delay\n");
	    // Based on ELL1Hmodel.C
	    
	    //h3 = psr[p].param[param_h3].val[0];
	    h3 = getParameterValue( &psr[p], param_h3, 0 );
	    
	    // Determine fw10 mode
	    if( psr[p].param[param_h4].paramSet[0] == 1 ){
	      //	      h4 = psr[p].param[param_h4].val[0];
	      h4 = getParameterValue( &psr[p], param_h4, 0 );
	      // mode 2 or 3 take preference over mode 1 as they are more stable
	      if( psr[p].param[param_nharm].paramSet[0] == 1 ){
		nharm = (int)psr[p].param[param_nharm].val[0];
		//nharm = (int)getParameterValue( &psr[p], param_nharm, 0 );
		if( nharm > 4 )
		  mode = 3;
		else
		  mode = 2;
	      }
	      if( psr[p].param[param_stig].paramSet[0] == 1 ){
		// Conflict. Unsure whether to select mode 1 or modes 2/3, so will default
		// to the most stable one.
		printf( "WARNING! You specified both H4 and STIG.\n" );
		printf( "We will ignore STIG and perform the approx. H4 fit instead.\n" );
		printf( "If you want to perform the exact fit for H3 and STIG, then " );
		printf( "please remove H4 from your parameter file.\n");
	      }
	      // Have H3, H4, but no NHARM
	      mode = 2;
	    }else{ 
	      // Have H3, but no H4
	      if( psr[p].param[param_stig].paramSet[0] == 1 ){
		//		stig = psr[p].param[param_stig].val[0];
		stig = getParameterValue( &psr[p], param_stig, 0 );
		mode = 1;
	      }else{
		mode = 0;
		h4 = 0;
		nharm = 3;
	      }
	    }// fw10 mode determined.
	    
	    // Define sin(i) and m2 for calculation of the orbital phases etc.
	    if( mode == 1 ){
	      // fw10, Eq. 22:
	      si = 2.0 * stig / ( 1.0 + pow( stig, 2.0 ) );
	      // fw10, Eq. 20:
	      m2 = h3 / pow( stig, 3.0 ); // Shapiro r, not just M2.
	      
	      if( si > 1.0 ){
		displayMsg(1,"BIN1",
			   "SIN I > 1.0, setting to 1: should probably use DDS model",
			   "",psr[p].noWarnings);
		si = 1.0;
		psr[p].param[param_sini].val[0] = 1.0L;
	      }
	    }else if( mode == 2 || mode == 3 ){
	      // fw10, Eq. 25:
	      si = 2.0 * h3 * h4 / ( h3 * h3 + h4 * h4 );
	      // fw10, Eq. 26:
	      m2 = pow( h3, 4.0 ) / pow( h4, 3.0 );
	      if( si > 1.0 ){
		displayMsg(1,"BIN1",
			   "SIN I > 1.0, setting to 1: should probably use DDS model",
			   "",psr[p].noWarnings);
		si = 1.0;
		psr[p].param[param_sini].val[0] = 1.0L;
	      }
	    }else if( mode == 0 ){
	      // Cannot determine m2 and/or sini. Will have to determine the
	      // Shapiro delay based on h3 alone.
	    }else{
	      printf( "This should not be possible. Go Away.\n" );
	      printf( "And tell someone about it: joris.verbiest@gmail.com, e.g.\n" );
	    }
	    brace=1-si*sin(phase);
	    dlogbr=log(brace);
	    
	    ecc = sqrt( eps1 * eps1 + eps2 * eps2 );
	    TrueAnom = phase;
	    //TrueAnom = 2.0 * atan2( sqrt( 1.0 + ecc ) * sin( phase / 2.0 ), 
	    //                       sqrt( 1.0 - ecc ) * cos( phase / 2.0 ) );
	    omega = atan2( eps1, eps2 );
	    //lgf = log( 1.0 + stig * stig - 2.0 * stig * sin( TrueAnom + omega ) );

	    fs = 1.0 + stig * stig - 2.0 * stig * sin( TrueAnom );
	    lgf = log( fs );
	    lsc = lgf + 2.0 * stig * sin( TrueAnom ) - stig * stig * cos( 2.0 * TrueAnom );
	    
	    if( mode == 0 ){
	      // mode 0: only h3 is known. 
	      ds = -4.0 / 3.0 * h3 * sin( 3.0 * TrueAnom );
	    }else if( mode == 1 ){
	      ds = -2.0 *m2* lsc;
	      //ds = -2.0 * m2 * dlogbr;
	    }else{ // modes 2 and 3
	      ds = calcDH( TrueAnom, h3, h4, nharm, 0 );
	    }
	    */
	  //	}   // End of Shapiro delay parameters set if
	  //	  else{
	    dlogbr=log(brace);
	    ds=-2*m2*dlogbr;        /* Equation 26 */
	    //	  }
	  
        }
      else
        {
          printf("ERROR T2model: Require eccentricity set or EPS1/EPS2 parameters for companion %d\n",com+1);
	  printf("ecc=%lf, eps1=%lf, eps2=%lf\n", ecc, eps1, eps2);
          exit(1);
        }
      // printf("T2: %g %g %g %g %g %g %g\n",brace,phase,a0,b0,dlogbr,ds,m2);

      /* Now compute d2bar, the orbital time correction in DD equation 42. */
      /* Equation 52 */
      if (onemecu != 0.0)
        {
          d2bar=dre*(1-anhat*drep+allTerms*pow(anhat,2)*
                     (pow(drep,2) + 0.5*dre*drepp - 
                      0.5*ecc*su*dre*drep/onemecu)) + allTerms*(ds+da+DAOP+DSR
                                                                + DOP);
        }
      else
        {
          d2bar=dre*(1-anhat*drep+allTerms*pow(anhat,2)*
                     (pow(drep,2) + 0.5*dre*drepp))
            + allTerms*(ds+da+DAOP+DSR+DOP);
        }    
      //      printf("T2a: %g %g %g %g %g drepp=%g ecc=%g su =%g ome =%g ds = %g da = %g %g %g\n",(double)d2bar,(double)dre,(double)anhat,(double)drep,(double)allTerms,(double)drepp,(double)ecc,(double)su,(double)onemecu,(double)ds,(double)da,(double)DAOP,(double)DSR);
      torb-=d2bar;                                  /* Equation 42  */
    }
  // Removed all the param != -1 statements
  return torb;
}

/* ************************************************************************
   m2 - solves mass function for m2, using the Newton-Raphson method

where:  mf = mass function
m1 = primary mass
si = sin(i) (i = inclination angle)

solves: (m1+m2)^2 = (m2*si)^3 / mf

returns -1 on error

WVS Jan 2000
 ************************************************************************ */

double solve_mass_function_m2_double(long double mf, long double sini, long double m1)
{
   double guess = m1;
   double dx = 0.0;
   double eq = 0.0;
   double deq_dm2 = 0.0;
   int gi = 0;

   for (gi=0; gi<10000; gi++) {
	  eq = pow(m1+guess,2) - pow(guess*sini,3) / mf;
	  deq_dm2 = 2.0*(m1+guess) - 3.0 * pow(guess*sini,2) / mf;

	  dx = eq / deq_dm2;
	  guess -= dx;

	  if (fabs (dx) <= fabs(guess)*1e-10)
		 return guess;
   }
   printerror(0, "ERROR solve_mass_function_m2: maximum iterations exceeded\n");
   return -1.0;
}

long double solve_mass_function_m2(long double mf, long double sini, long double m1)
{
   long double guess = m1;
   long double dx = 0.0;
   long double eq = 0.0;
   long double deq_dm2 = 0.0;
   long double guesssini, guesssini2;
   int gi = 0;

   for (gi=0; gi<10000; gi++) {
     guesssini = guess*sini;
     guesssini2 = guesssini*guesssini;
     eq = (m1+guess)*(m1+guess) - guesssini2*guesssini / mf;
     deq_dm2 = 2.0*(m1+guess) - 3.0 * guesssini2 / mf;

     dx = eq / deq_dm2;
     guess -= dx;

     if (fabsl (dx) <= fabsl(guess)*1e-10)
       return guess;
   }
   printerror(0, "ERROR solve_mass_function_m2: maximum iterations exceeded\n");
   return -1.0;
}

/* Compute eccentric anomaly u by iterating Kepler's equation */	 
/*  Compute eccentric anomaly u by iterating Kepler's equation if
    eccentricity is set. The equation is solved using a Newton-Raphson
    technique and the S9 * starting value in Odell & Gooding 1986
    CeMec 38 307 */
// Patrick: changed du in a long double to avoid the following iteration getting stuck, also defined a new u variable to use temporarily. Also added counter to avoid getting stuck.
// Patrick: Added hyperbolic version
// Patrick: Parabolic orbits: mean anomaly gives true anomaly directly rather than an "eccentric anomaly"
// phase = (2pi/PB)*((t-t0) - 0.5*(pbdot/pb)*(t-t0)**2)
long double getEccentricAnomaly(long double phase, long double ecc)
{
  long double du, new_u;
  long nrsteps;
  nrsteps = 0;

  if(ecc < 0) {
    printerror(0, "ERROR: In getEccentricAnomaly() Negative eccentricities are not possible");
    exit(0);
  }else if(ecc < 1.0) {
    /*	  u=phase+ecc*sin(phase)*(1.0+ecc*cos(phase));*/
    new_u = phase+ecc*sin(phase)/sqrt(1.0-2.0*ecc*cos(phase)+ecc*ecc);
    //	  fprintf(stderr, "Starting Kepler's equation\n");
    do {
      du=(phase-(new_u-ecc*sinl(new_u)))/(1.0-ecc*cosl(new_u));
      new_u += du;
      nrsteps++;
      if(nrsteps == 1000) {
	printwarning(0, "WARNING: In getEccentricAnomaly() iteration of Kepler's equation didn't reached required precision (e=%Lf)", ecc);
	break;
      }
    }while (fabsl(du)>1.0e-14);
  }else if(ecc > 1.0) {
    //    new_u = powl(6.0*phase, 1.0/3.0);  // Not sure what a suitable start condition is. (6*phase)^(1/3) has been suggested?
    new_u = phase;   // Prussing 1977 - J. Astronaut. Sci., 25, 123 - See also Toshio Fukushima - A method solving Kepler's equation for hyperbolic case - Celestial Mechanics and Dynamical Astronomy - 1997, 68, 121
    do {
      du=(phase-(-new_u+ecc*sinhl(new_u)))/(ecc*coshl(new_u)-1.0);
      new_u += du;
      nrsteps++;
      if(nrsteps == 1000) {
	printwarning(0, "WARNING: In getEccentricAnomaly() iteration of Kepler's equation didn't reached required precision (e=%Lf)", ecc);
	break;
      }
    }while (fabsl(du)>1.0e-14);
  }else {  
    //  }else if(ecc == 1.0) { // Should always be ecc==1.0 at this point, but explicit if results in a gcc warning
    new_u = 3.0*phase;
    new_u = new_u + sqrtl(new_u*new_u+1.0);
    new_u = powl(new_u, 1.0/3.0)-powl(new_u, -1.0/3.0);
    new_u = 2.0*atanl(new_u);
  }
  //	  fprintf(stderr, "Done with Kepler's equation\n");
  return new_u;
}


