#include <math.h>

#include <psrsalsa.h>

#define SecPerDay           86400.0         /* 24.0*3600.0 */

// This is the same as T2model()


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
//
// Removed post-Keplerian parameters, but added hyperbolic orbits and parabolic orbits.
//
// select_companion = -1:  add effects for all companions
// select_companion != -1: Only consider selected companion (counting from zero)
// calc_deriv = 0: return time correction
// calc_deriv = 1: return derivative of time correction to a
// calc_deriv = 2: return derivative of time correction to Pb
// calc_deriv = 3: return derivative of time correction to T0
// calc_deriv = 4: return derivative of time correction to omega
// calc_deriv = 5: return derivative of time correction to eps
double T2model_hyper(long double bbat, int nCompanion, int calc_deriv, int select_companion, binary_def *companion)
{
  double omegaB;
  double pb, pbdot;
  double rad2deg = 180.0/M_PI;
  double tt0,t0,x,xdot, asini, ecc;
  double phase,u;
  double orbits;
  int    norbits;
  double omega;
  double sinE, cosE, ecosE_min_1, alpha, beta, delta, delta1, delta2, nhat;
  double d2bar,torb;
  int    com;
  int    allTerms=1;            /* = 0 to emulate BT model */   // Assume we can approximate things in the BT model. I think it fails if [2pi/pb(in sec)]/(1-ecc*cos(E)) is not << 1, hence it could go wrong if ecc ~ 1 and E~0=near periastron

  torb = 0.0;

  for (com = 0; com < nCompanion; com++) {      
    if(select_companion >= 0)
      com = select_companion;
    pb  = companion[com].pb*SecPerDay;
    t0  = companion[com].t0;
    ecc = companion[com].ecc; 
    //    printf("XXXXX: calc_deriv=%d ecc=%e\n", calc_deriv, ecc);
    //    fflush(stdout);
    omega = companion[com].om/rad2deg;
    x   = companion[com].a1;   // Note: for parabolic orbits this really is interpreted to be "semi-latus rectum" = a(1-e^2). It is therefore not smoothly connecting to hyperbolic/elliptic orbits. Could abandon a for p in all models, which should be a more sensible approach. However, for compatibility this has not been done.
    xdot  = companion[com].a1dot; 
    pbdot = companion[com].pbdot;

    omegaB    = 2.0*M_PI/pb;   //This is the mean angular frequency

    //tt0 = TOA - reference epoch
    if (t0 != 0) {    
      tt0 = (bbat-t0)*SecPerDay;    
    }else {
      printf("ERROR T2model_hyper: T0 needs to be set in the parameter file\n");
      exit(1);
    }
    
    asini  = x + xdot*tt0;


    /* Do some checks */
    if (ecc < 0.0) {
      printerror(0, "ERROR T2model_hyper: problem with eccentricity = %e, quiting.", ecc);
      exit(0);
      //          psr[p].param[param_ecc].val[com]=0.0;
      //          ecc = 0.0;
    }
    if(calc_deriv > 0) {
      if(xdot != 0.0) {
	printerror(0, "ERROR T2model_hyper: Derivatives are not tested with an adot, quiting.");
	exit(0);
      }
      if(pbdot != 0.0) {
	printerror(0, "ERROR T2model_hyper: Derivatives are not tested with a pbdot, quiting.");
      }
    }
    
    /* Obtain number of orbits in tt0 */
    if(ecc < 1) {
      orbits  = tt0/pb - 0.5*(pbdot)*pow(tt0/pb,2);
      norbits = (int)orbits;
      if (orbits<0.0) norbits--;
      /* Obtain phase of orbit */
      phase=2.0*M_PI*(orbits-norbits);
    }else {
      if(pbdot != 0) {
	printf("ERROR T2model_hyper: PBDOT undefined for hyperbolic/parabolic orbits\n");
	exit(1);	
      }
      phase=2.0*M_PI*tt0/pb;  // Mean anomaly doesn't need orbit subtraction as there is only one orbit
    }

    //    printf("XXXXXX HYPER??? ecc=%lf\n", ecc);
    if(companion[com].ecc_defined) {
      /* Compute eccentric anomaly u by iterating Kepler's equation */
      // But for parabolic orbits the true anomaly is returned, which is determined analytically
      u = getEccentricAnomaly(phase, ecc);
      if(ecc < 1.0) {  // Elliptical orbit
	sinE = sin(u);
	cosE = cos(u);
	ecosE_min_1 = 1.0-ecc*cosE;
	alpha = asini*sin(omega);                                   /* Equation 46  */    // = a1*sin omega
	beta = asini*sqrt(1.0-ecc*ecc)*cos(omega);                 /* Equation 47  */    // = a1*sqrt(1-ecc**2)*cos omega
	delta=alpha*(cosE-ecc) + beta*sinE;                    /* Equation 48  */    // = This is the "A" term from my notes
	delta1=-alpha*sinE + beta*cosE;                       /* Equation 49  */    
	delta2=-alpha*cosE - beta*sinE;                      /* Equation 50  */    // 
	nhat=omegaB/ecosE_min_1;                             /* Equation 51  */    // = MeanAnomaly/(1-ecc*cos(E)
	/* Now compute d2bar, the orbital time correction in DD equation 42. */
	/* Equation 52 */
	if (ecosE_min_1 != 0.0) {  // Only a problem if ecc=1, which should be handled separately anyways.
	  if(calc_deriv == 0) {
	    d2bar=delta*(1-nhat*delta1+allTerms*nhat*nhat*(delta1*delta1 + 0.5*delta*delta2 - 0.5*ecc*sinE*delta*delta1/ecosE_min_1));
	  }else if(calc_deriv == 1) {
	    d2bar=(delta/x)*(1-2.0*nhat*delta1+allTerms*3.0*nhat*nhat*(delta1*delta1 + 0.5*delta*delta2 - 0.5*ecc*sinE*delta*delta1/ecosE_min_1));
	  }else if(calc_deriv == 2) {   // Derivative to Pb
	    double n2 = nhat*nhat;
	    double ratio = -ecc*sinE/ecosE_min_1;
	    double ratio2 = -ecc*cosE/ecosE_min_1;
	    // Calculate derivative to E
	    d2bar = delta1*(1.0 - nhat*delta1 + n2*(delta1*delta1+3.0*delta*(delta2-delta1*ratio)));
	    d2bar += delta*(nhat*ratio*(delta1-1.5*nhat*delta*delta2)-nhat*delta2
			    -0.5*n2*delta*delta1*(1.0+ratio2-3.0*ratio*ratio));
	    // Multiply with derivative to M
	    d2bar /= -ecosE_min_1;
	    // Multiply with derivative to Pb
	    d2bar *= 2.0*M_PI*tt0*SecPerDay/(pb*pb);
	  }else if(calc_deriv == 3) {   // Derivative to T0
	    double n2 = nhat*nhat;
	    double ratio = -ecc*sinE/ecosE_min_1;
	    double ratio2 = -ecc*cosE/ecosE_min_1;
	    // Calculate derivative to E
	    d2bar = delta1*(1.0 - nhat*delta1 + n2*(delta1*delta1+3.0*delta*(delta2-delta1*ratio)));
	    d2bar += delta*(nhat*ratio*(delta1-1.5*nhat*delta*delta2)-nhat*delta2
			    -0.5*n2*delta*delta1*(1.0+ratio2-3.0*ratio*ratio));
	    // Multiply with derivative to M  (so far the same as for the derivative to Pb)
	    d2bar /= -ecosE_min_1;
	    // Multiply with derivative to T0
	    d2bar *= omegaB*SecPerDay;
	  }else if(calc_deriv == 4) {   // Derivative to omega
	    double tmp = sqrt(1.0-ecc*ecc);
	    double betaprime = beta/tmp;
	    double alphaprime = alpha*tmp;
	    double n2 = nhat*nhat;
	    double ratio = -ecc*sinE/ecosE_min_1;
	    double tmp2;
	    // Calculate derivative to F (same as for Pb)
	    d2bar = (betaprime*(cosE-ecc)-alphaprime*sinE);
	    //	  printf("XXXX a=%le b=%le\n", alpha, beta);
	    //	  printf("XXXX %le %le %le %le %le\n", betaprime, ecc, cosE, alphaprime, sinE);
	    //	  printf("XXXX %le\n", d2bar);
	    d2bar *= (1.0 - nhat*delta1 + n2*(delta1*delta1+0.5*delta*(delta2-delta1*ratio)));
	    tmp = (-betaprime*sinE-alphaprime*cosE);
	    tmp2 = (betaprime*(cosE-ecc)-alphaprime*sinE);
	    d2bar += delta*(nhat*(betaprime*sinE+alphaprime*cosE)+2.0*n2*delta1*tmp
			    +0.5*n2*(tmp2*delta2 + delta*(-betaprime*cosE+alphaprime*sinE)
				     - (tmp2*delta1 + delta*tmp*ratio)));
	    d2bar /= rad2deg;
	  }else if(calc_deriv == 5) {   // Derivative to ecc
	    double betaprime = -beta*ecc/(1.0-ecc*ecc);
	    double n2 = nhat*nhat;
	    double Eprime = -sinE/ecosE_min_1;
	    double eEprime = ecc*Eprime;
	    double tmp;
	    // Calculate derivative to F (same as for Pb)
	    d2bar = -(-alpha+Eprime*delta1+betaprime*sinE);
	    d2bar *= (1.0 - nhat*delta1 + n2*(delta1*delta1+0.5*delta*(delta2-delta1*eEprime)));
	    tmp = (eEprime*sinE-cosE)/ecosE_min_1;
	    d2bar -= delta*( -nhat*(Eprime*delta2+betaprime*cosE) - nhat*tmp*delta1
			     +2.0*n2*tmp*(delta1*delta1+0.5*delta*(delta2-delta1*eEprime))
			     +n2*(2.0*delta1*(Eprime*delta2+betaprime*cosE)
				  +0.5*(-alpha+Eprime*delta1+betaprime*sinE)*(delta2-delta1*eEprime)
				  +0.5*delta*(-Eprime*delta1-betaprime*sinE-(Eprime*delta2+betaprime*cosE)*eEprime
					      -delta1*Eprime*(1.0-2.0*ecc*cosE/ecosE_min_1 - eEprime*eEprime)
					      )));
	    d2bar *= -1.0;
	  }else {
	    printf("ERROR T2model_hyper: Requested derivative is not supported (calc_deriv = %d for elliptical orbit).\n", calc_deriv);
	    exit(1);
	  }
        }else {
	  if(calc_deriv == 0) {
	    d2bar=delta*(1-nhat*delta1+allTerms*nhat*nhat*(delta1*delta1 + 0.5*delta*delta2));
	  }else if(calc_deriv == 1) {
	    d2bar=(delta/x)*(1-2.0*nhat*delta1+allTerms*3.0*nhat*nhat*(delta1*delta1 + 0.5*delta*delta2));
	  }else {
	    printf("ERROR T2model_hyper: Requested derivative is not supported (calc_deriv = %d for elliptical orbit with e cos E = 1).\n", calc_deriv);
	    exit(1);
	  }
	}          
      }else if(ecc == 1.0) {   // Up to first order, unlike hyperbolic/excentric case
	double sinu, cosu, E, tan_half_u;
	sinu = sin(u);   // u is true anomaly, phase=mean anomaly
	cosu = cos(u);
	E = 3.0*phase+sqrt(9.0*phase*phase+1.0);
	alpha = asini/(1.0+cosu);
	beta = alpha*cos(omega);
	alpha *= sin(omega);
	delta = alpha*cosu + beta*sinu;                    // Equation 48      // = This is the "A" term from my notes
	tan_half_u = tan(0.5*u);
	if(calc_deriv == 0) {
	  d2bar = delta*(1.0-(beta*cosu-alpha*sinu+delta*sinu/(1.0+cosu))*2.0*omegaB*(pow(E, 1.0/3.0)+pow(E, -1.0/3.0))/((1.0+tan_half_u*tan_half_u)*(E-3.0*phase)));
	}else {
	  printf("ERROR T2model_hyper: Requested derivative is not supported (calc_deriv = %d for parabolic orbit).\n", calc_deriv);
	  exit(1);
	}
	/*
	sinE = sin(u);
	cosE = cos(u);
	alpha = asini*sin(omega);                                   // Equation 46      // = a1*sin omega
	beta = asini*cos(omega);                 // Equation 47      // = a1*cos omega
	delta = alpha*(cosE) + beta*sinE;                    // Equation 48      // = This is the "A" term from my notes
	delta /= 1.0+cosE;
	if(calc_deriv == 0) {
	  d2bar = delta;
	}else if(calc_deriv == 1) {
	  d2bar = delta/x;
	}else {
	  printf("ERROR T2model_hyper: Requested derivative is not supported (calc_deriv = %d for parabolic orbit).\n", calc_deriv);
	  exit(1);
	}
*/
	//	printf("XXXXXX HYPER!!!!\n");
      }else {
	sinE = sinh(u);
	cosE = cosh(u);
	ecosE_min_1 = ecc*cosE-1.0;       // Denominator of nhat
	alpha = asini*sin(omega);                                   /* Equation 46  */    // = a1*sin omega
	beta = asini*sqrt(ecc*ecc-1.0)*cos(omega);                 /* Equation 47  */    // = a1*sqrt(1-ecc**2)*cos omega
	delta = alpha*(ecc-cosE) + beta*sinE;                    /* Equation 48  */    // = This is the "A" term from my notes
	delta1 = -alpha*sinE + beta*cosE;                       /* Equation 49  */    
	delta2 = -alpha*cosE + beta*sinE;                      /* Equation 50  */    // 
	nhat = omegaB/ecosE_min_1;                             /* Equation 51  */    // = MeanAnomaly/(1-ecc*cos(E)
	if(calc_deriv == 0) {
	  d2bar = delta*(1.0-nhat*delta1+allTerms*nhat*nhat*(delta1*delta1 - 0.5*delta*delta2 - 0.5*ecc*sinE*delta*delta1/ecosE_min_1));
	}else if(calc_deriv == 1) {   // Derivative to a
	  //	  printf("XXXXX com=%d, mjd=%Le, phase=%e, tt0=%e, pb=%e, x=%e, u=%e, delta=%e, alpha=%e, beta=%e, ecc=%e, cosE=%e\n", com, bbat, phase, tt0, pb, x, u, delta, alpha, beta, ecc, cosE);
	  d2bar = (delta/x)*(1.0-2.0*nhat*delta1+allTerms*3.0*nhat*nhat*(delta1*delta1 - 0.5*delta*delta2 - 0.5*ecc*sinE*delta*delta1/ecosE_min_1));
	}else if(calc_deriv == 2) {   // Derivative to Pb
	  double n2 = nhat*nhat;
	  double ratio = ecc*sinE/ecosE_min_1;
	  double ratio2 = ecc*cosE/ecosE_min_1;
	  // Calculate derivative to F
	  d2bar = delta1*(1.0 - nhat*delta1 + n2*(delta1*delta1+delta*(delta2-3.0*delta1*ratio)));
	  d2bar += delta*(nhat*(ratio*delta1+delta2)
			  -0.5*n2*delta*(-delta2*ratio+delta1*(1.0+ratio2-3.0*ratio*ratio)));
	  // Multiply with derivative to M
	  d2bar /= ecosE_min_1;
	  // Multiply with derivative to Pb
	  d2bar *= -2.0*M_PI*tt0*SecPerDay/(pb*pb);
	}else if(calc_deriv == 3) {   // Derivative to T0
	  double n2 = nhat*nhat;
	  double ratio = ecc*sinE/ecosE_min_1;
	  double ratio2 = ecc*cosE/ecosE_min_1;
	  // Calculate derivative to F (same as for Pb)
	  d2bar = delta1*(1.0 - nhat*delta1 + n2*(delta1*delta1+delta*(delta2-3.0*delta1*ratio)));
	  d2bar += delta*(nhat*(ratio*delta1+delta2)
			  -0.5*n2*delta*(-delta2*ratio+delta1*(1.0+ratio2-3.0*ratio*ratio)));
	  // Multiply with derivative to M (same as for Pb)
	  d2bar /= ecosE_min_1;
	  // Multiply with derivative to T0
	  d2bar *= -omegaB*SecPerDay;
	}else if(calc_deriv == 4) {   // Derivative to omega
	  double tmp = sqrt(ecc*ecc-1.0);
	  double betaprime = beta/tmp;
	  double alphaprime = alpha*tmp;
	  double n2 = nhat*nhat;
	  double ratio = ecc*sinE/ecosE_min_1;
	  double tmp2;
	  // Calculate derivative to F (same as for Pb)
	  d2bar = (betaprime*(ecc-cosE)-alphaprime*sinE);
	  //	  printf("XXXX a=%le b=%le\n", alpha, beta);
	  //	  printf("XXXX %le %le %le %le %le\n", betaprime, ecc, cosE, alphaprime, sinE);
	  //	  printf("XXXX %le\n", d2bar);
	  d2bar *= (1.0 - nhat*delta1 + n2*(delta1*delta1-0.5*delta*(delta2+delta1*ratio)));
	  tmp = (-betaprime*sinE-alphaprime*cosE);
	  tmp2 = (betaprime*(ecc-cosE)-alphaprime*sinE);
	  d2bar += delta*(nhat*(betaprime*sinE+alphaprime*cosE)+2.0*n2*delta1*tmp
			  -0.5*n2*(tmp2*delta2 + delta*(-betaprime*cosE-alphaprime*sinE)
				   + tmp2*delta1 + delta*tmp*ratio));
	  d2bar /= rad2deg;
	}else if(calc_deriv == 5) {   // Derivative to ecc
	  double betaprime = beta*ecc/(ecc*ecc-1.0);
	  double n2 = nhat*nhat;
	  double Fprime = -sinE/ecosE_min_1;
	  double eFprime = ecc*Fprime;
	  double tmp;
	  // Calculate derivative to F (same as for Pb)
	  d2bar = -(alpha+Fprime*delta1+betaprime*sinE);
	  d2bar *= (1.0 - nhat*delta1 + n2*(delta1*delta1-0.5*delta*(delta2+delta1*eFprime)));
	  tmp = (eFprime*sinE+cosE)/ecosE_min_1;
	  d2bar -= delta*( -nhat*(Fprime*delta2+betaprime*cosE) + nhat*tmp*delta1
			   -2.0*n2*tmp*(delta1*delta1-0.5*delta*(delta2+delta1*eFprime))
			   +n2*(2.0*delta1*(Fprime*delta2+betaprime*cosE)
				-0.5*(alpha+Fprime*delta1+betaprime*sinE)*(delta2+delta1*eFprime)
				-0.5*delta*(Fprime*delta1+betaprime*sinE+(Fprime*delta2+betaprime*cosE)*eFprime
					    //					    +delta1*((sinE+eFprime*cosE)/ecosE_min_1 - eFprime*(cosE+eFprime*sinE)/(ecosE_min_1*ecosE_min_1))
					    -delta1*Fprime*(1.0-2.0*ecc*cosE/ecosE_min_1 - eFprime*eFprime)
					    )));
	  d2bar *= -1.0;
	}else {
	  printf("ERROR T2model_hyper: Requested derivative is not supported.\n");
	  exit(1);
	}
      }
      
      
    }else {
      printf("ERROR T2model_hyper: Require eccentricity to be set for companion %d (EPS1/EPS2 parameters not supported)\n",com+1);
      exit(1);
    }

    torb-=d2bar;                                  /* Equation 42  */
    if(select_companion >= 0)
      break;
  }
  return torb;
}



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

//Differences with BTmodel:
//  - Both hyperbolic and elliptical orbits are allowed
//  - Switched to generalised Kepler's equation solver
//  - No post-Keplerian parameters are handled
//  - Multiple companions are supported
double BTmodel_hyper(long double bbat, int nCompanion, binary_def *companion)
{
  double torb = 0;
  double tt0;
  double orbits;
  long   norbits, compnr;
  double phase;
  double bige,tt,som,com;
  double alpha,beta,sbe,cbe,q,r,s;


  for (compnr = 0; compnr < nCompanion; compnr++) {      

    long double pb    = companion[compnr].pb*86400.0;
    long double ecc   = companion[compnr].ecc;
    long double omega = companion[compnr].om*M_PI/180.0;   // omega (longitude of periastron) in radians
    
    tt0 = (bbat - companion[compnr].t0)*86400.0;  // time since epoch in sec
    
    if (ecc < 0.0 || ecc == 1.0) {
      printerror(0, "ERROR BTmodel_hyper: problem with eccentricity = %Lg\n", ecc);
      exit(0);
    }
    
    /* Should ct be the barycentric arrival time? -- check bnrybt.f */
    if(ecc < 1) {
      orbits = tt0/pb;   // time since epoch in orbits
      norbits = (int)orbits;
      if (orbits < 0.0) norbits--;
      phase = 2.0*M_PI * (orbits-norbits);
    }else {
      phase = 2.0*M_PI*tt0/pb;  // Mean anomaly doesn't need orbit subtraction
    }
    
    bige = getEccentricAnomaly(phase, ecc);

    if(ecc < 1)
      tt = 1.0-ecc*ecc;
    else
      tt = ecc*ecc-1.0;
    som = sin(omega);
    com = cos(omega);
    
    alpha = companion[compnr].a1*som;           // Eq. 2.31
    beta = companion[compnr].a1*com*sqrt(tt);   // Eq. 2.31
    if(ecc < 1) {
      sbe = sin(bige);
      cbe = cos(bige);
      q = alpha * (cbe-ecc) + beta*sbe;   // - first two corrections in 2.33
      r = -alpha*sbe + beta*cbe;         // - first term in numerator of last term in 2.33
      s = 1.0/(1.0-ecc*cbe);             // denominator
    }else {
      sbe = sinh(bige);
      cbe = cosh(bige);
      q = alpha * (ecc-cbe) + beta*sbe;   // - first two corrections in 2.33
      r = -alpha*sbe + beta*cbe;         // - first term in numerator of last term in 2.33
      s = 1.0/(ecc*cbe-1.0);             // denominator
    }
    
    torb += -q+(2*M_PI/pb)*q*r*s;      // This is Eq. 2.33
  }

  return torb;

}

