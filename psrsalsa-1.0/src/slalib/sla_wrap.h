//START REGION RELEASE

/*
SLA_ALTAZ - Velocities etc. for Altazimuth Mount   

ACTION: Positions, velocities and accelerations for an altazimuth telescope mount tracking a star (double precision).
CALL: CALL sla_ALTAZ ( HA, DEC, PHI, AZ, AZD, AZDD, EL, ELD, ELDD, PA, PAD, PADD)
GIVEN:
    HA 	D 	hour angle
    DEC 	D 	declination
    PHI 	D 	observatory latitude
RETURNED:
    AZ 	D 	azimuth
    AZD 	D 	azimuth velocity
    AZDD 	D 	azimuth acceleration
    EL 	D 	elevation
    ELD 	D 	elevation velocity
    ELDD 	D 	elevation acceleration
    PA 	D 	parallactic angle
    PAD 	D 	parallactic angle velocity
    PADD 	D 	parallactic angle acceleration

NOTES:
    1. Natural units are used throughout. HA, DEC, PHI, AZ, EL and ZD are in radians. The velocities and accelerations assume constant declination and constant rate of change of hour angle (as for tracking a star); the units of AZD, ELD and PAD are radians per radian of HA, while the units of AZDD, ELDD and PADD are radians per radian of HA squared. To convert into practical degree- and second-based units:

        angles 	       *360/2pi 	       ->	degrees
        velocities     *(2pi/86400)*(360/2pi)  -> 	degree/sec
        accelerations *(2pi/86400)^2*(360/2pi) -> 	degree/sec/sec

        Note that the seconds here are sidereal rather than SI. One sidereal second is about 0.99727 SI seconds.

        The velocity and acceleration factors assume the sidereal tracking case. Their respective numerical values are (exactly) 1/240 and (approximately) 1/3300236.9. 
    2. Azimuth is returned in the range [0,2pi]; north is zero, and east is +pi/2. Elevation and parallactic angle are returned in the range +-pi/2. Position angle is +ve for a star west of the meridian and is the angle NP-star-zenith. 
    3. The latitude is geodetic as opposed to geocentric. The hour angle and declination are topocentric. Refraction and deficiencies in the telescope mounting are ignored. The purpose of the routine is to give the general form of the quantities. The details of a real telescope could profoundly change the results, especially close to the zenith. 
    4. No range checking of arguments is carried out. 
    5. In applications which involve many such calculations, rather than calling the present routine it will be more efficient to use inline code, having previously computed fixed terms such as sine and cosine of latitude, and (for tracking a star) sine and cosine of declination. 
 */
void csla_altaz(double ha, double dec, double phi, double *az, double *azd, double *azdd, double *el, double *eld, double *eldd, double *pa, double *pad, double *padd);

//START REGION DEVELOP

/*
SLA_CLDJ - Calendar to MJD   

ACTION: Gregorian Calendar to Modified Julian Date. 
CALL:   CALL sla_CLDJ (IY, IM, ID, DJM, J)
GIVEN:  IY,IM,ID 	I 	year, month, day in Gregorian calendar
RETURNED:
    DJM 	D 	modified Julian Date (JD-2400000.5) for 0h
    J 	I 	status:
    		0 = OK
    		1 = bad year
    		2 = bad month
    		3 = bad day
NOTES:
    1. When an invalid year or month is supplied (status J = 1 or 2) the MJD is not computed. When an invalid day is supplied (status J = 3) the MJD is computed. 
    2. The year must be -4699 (i.e. 4700BC) or later. For year nBC use IY = -(n-1). 
    3. An alternative to the present routine is sla_CALDJ, which accepts a year with the century missing. 

REFERENCE:
    The algorithm is derived from that of Hatcher, Q.Jl.R.astr.Soc. (1984) 25, 53-55. 
*/
void csla_cldj(int iy, int im, int id, double *djm, int *j);


//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_DCC2S - Cartesian to Spherical   

ACTION: Cartesian coordinates to spherical coordinates (double precision). 
CALL:   CALL sla_DCC2S (V, A, B)
GIVEN:
    V 	D(3) 	x,y,z vector
RETURNED:
    A,B 	D 	spherical coordinates in radians
NOTES:

    1. The spherical coordinates are longitude (+ve anticlockwise looking from the +ve latitude pole) and latitude. The Cartesian coordinates are right handed, with the x-axis at zero longitude and latitude, and the z-axis at the +ve latitude pole. 
    2. If V is null, zero A and B are returned. 
    3. At either pole, zero A is returned. 
 */
void csla_dcc2s(double v[3], double *a, double *b);

//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_DCS2C - Spherical to Cartesian   

ACTION: Spherical coordinates to Cartesian coordinates (double precision). 
CALL:   CALL sla_DCS2C (A, B, V)
GIVEN:
    A,B 	D 	spherical coordinates in radians: alpha,delta etc.
RETURNED:
    V 	D(3) 	x,y,z unit vector

NOTE: The spherical coordinates are longitude (+ve anticlockwise looking from the +ve latitude pole) and latitude. The Cartesian coordinates are right handed, with the x-axis at zero longitude and latitude, and the z-axis at the +ve latitude pole. 
*/
void csla_dcs2c(double a, double b, double v[3]);

//START REGION DEVELOP

/*
SLA_DE2H - h,delta to Az,El

ACTION: Equatorial to horizon coordinates (double precision).
CALL:   CALL sla_DE2H (HA, DEC, PHI, AZ, EL)
GIVEN:
    HA 	D 	hour angle (radians)
    DEC 	D 	declination (radians)
    PHI 	D 	latitude (radians)
RETURNED:
    AZ 	D 	azimuth (radians)
    EL 	D 	elevation (radians)

NOTES:
    1. Azimuth is returned in the range 0-2pi; north is zero, and east is +pi/2. Elevation is returned in the range +-pi. 
    2. The latitude must be geodetic. In critical applications, corrections for polar motion should be applied. 
    3. In some applications it will be important to specify the correct type of hour angle and declination in order to produce the required type of azimuth and elevation. In particular, it may be important to distinguish between elevation as affected by refraction, which would require the observed h,delta, and the elevation in vacuo, which would require the topocentric h,delta. If the effects of diurnal aberration can be neglected, the apparent h,delta may be used instead of the topocentric h,delta. 
    4. No range checking of arguments is carried out. 
    5. In applications which involve many such calculations, rather than calling the present routine it will be more efficient to use inline code, having previously computed fixed terms such as sine and cosine of latitude, and (for tracking a star) sine and cosine of declination. 
*/
void csla_de2h(double ha, double dec, double phi, double *az, double *el);

//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_DMXM - Multiply 3x3 Matrices

ACTION:    Product of two 3x3 matrices (double precision). 
CALL:      CALL sla_DMXM (A, B, C)
GIVEN:
    A 	D(3,3) 	matrix A
    B 	D(3,3) 	matrix B
RETURNED:
    C 	D(3,3) 	matrix result: AxB

NOTE: To comply with the ANSI Fortran 77 standard, A, B and C must be different arrays. The routine is, in fact, coded so as to work properly on the VAX and many other systems even if this rule is violated, something that is not, however, recommended. 
 */
void csla_dmxv(double dm[3][3], double va[3], double vb[3] );

//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_DRANRM - Put Angle into Range 0-2pi

ACTION: Normalize an angle into the range 0-2pi (double precision).

CALL: D = sla_DRANRM (ANGLE)
GIVEN:
    ANGLE 	D 	angle in radians
RETURNED:
    sla_DRANRM 	D 	ANGLE expressed in the range 0-2pi
 */
double csla_dranrm(double angle);

//START REGION DEVELOP

/*
SLA_DSEP - Angle Between 2 Points on Sphere   

ACTION:  Angle between two points on a sphere (double precision). 
CALL:    D = sla_DSEP (A1, B1, A2, B2)
GIVEN:
    A1,B1 	D 	spherical coordinates of one point (radians)
    A2,B2 	D 	spherical coordinates of the other point (radians)
RETURNED:
    sla_DSEP 	D 	angle between [A1,B1] and [A2,B2] in radians

NOTES:
    1. The spherical coordinates are right ascension and declination, longitude and latitude, etc., in radians. 
    2. The result is always positive. 
*/
double csla_dsep(double a1, double b1, double a2, double b2 );

/*
SLA_DTT - TT minus UTC   

ACTION: Compute Delta TT, the increment to be applied to Coordinated Universal Time UTC to give Terrestrial Time TT.
CALL:   D = sla_DTT (DJU)
GIVEN:
    DJU 	D 	UTC date as a modified JD (JD-2400000.5)
RETURNED:
    sla_DTT 	D 	TT-UTC in seconds

NOTES:
    1. The UTC is specified to be a date rather than a time to indicate that care needs to be taken not to specify an instant which lies within a leap second. Though in most cases UTC can include the fractional part, correct behaviour on the day of a leap second can be guaranteed only up to the end of the second 23h59m59s. 
    2. Pre 1972 January 1 a fixed value of 10 + ET-TAI is returned. 
    3. TT is one interpretation of the defunct timescale Ephemeris Time, ET. 
    4. See also the routine sla_DT, which roughly estimates ET-UT for historical epochs. 
*/
double csla_dtt(double dju );

//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_DVDV - Scalar Product   

ACTION:    Scalar product of two 3-vectors (double precision). 
CALL:      D = sla_DVDV (VA, VB)
GIVEN:
    VA 	D(3) 	first vector
    VB 	D(3) 	second vector
RETURNED:
    sla_DVDV 	D 	scalar product VA.VB
 */
double csla_dvdv(double va[3], double vb[3]);

//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_DVN - Normalize Vector   

ACTION:    Normalize a 3-vector, also giving the modulus (double precision). 
CALL:      CALL sla_DVN (V, UV, VM)
GIVEN:
    V 	D(3) 	vector
RETURNED:
    UV 	D(3) 	unit vector in direction of V
    VM 	D 	modulus of V

NOTE:
    If the modulus of V is zero, UV is set to zero as well. 
*/
double csla_dvn(double v[3], double uv[3], double *vm);

//START REGION DEVELOP

/*
SLA_ECLEQ - Ecliptic to Equatorial   

ACTION:    Transformation from ecliptic longitude and latitude to J2000.0 alpha,delta
CALL:      CALL sla_ECLEQ (DL, DB, DATE, DR, DD)
GIVEN:
    DL,DB 	D 	ecliptic longitude and latitude (mean of date, IAU 1980 theory, radians)
    DATE 	D 	TDB (formerly ET) as Modified Julian Date (JD-2400000.5)

RETURNED:
    DR,DD 	D 	J2000.0 mean alpha,delta (radians)
 */
void csla_ecleq(double elong, double elat, double mjd,
                double *raj, double *decj );


//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_EPB - MJD to Besselian Epoch   

ACTION: Conversion of Modified Julian Date to Besselian Epoch. 
CALL:   D = sla_EPB (DATE)
GIVEN:  DATE 	D 	Modified Julian Date (JD-2400000.5)
RETURNED: sla_EPB 	D 	Besselian Epoch

REFERENCE: Lieske, J.H., 1979, Astr.Astrophys. 73, 282. 
*/
double csla_epb(double date);

//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_EPJ - MJD to Julian Epoch   

ACTION:   Convert Modified Julian Date to Julian Epoch. 
CALL:     D = sla_EPJ (DATE)
GIVEN:    DATE 	D 	Modified Julian Date (JD-2400000.5)
RETURNED: sla_EPJ 	D 	Julian Epoch

REFERENCE: Lieske, J.H., 1979. Astr.Astrophys., 73, 282. 
*/
double csla_epj(double date );

//START REGION DEVELOP

/*
SLA_EQECL - J2000 alpha,delta to Ecliptic

ACTION:    Transformation from J2000.0 equatorial coordinates to ecliptic longitude and latitude.
CALL:     CALL sla_EQECL (DR, DD, DATE, DL, DB)
GIVEN:
    DR,DD 	D 	J2000.0 mean alpha,delta (radians)
    DATE 	D 	TDB (formerly ET) as Modified Julian Date (JD-2400000.5)
RETURNED:
    DL,DB 	D 	ecliptic longitude and latitude (mean of date, IAU 1980 theory, radians)
 */
void csla_eqecl(double raj, double decj, double mjd,
                double *elong, double *elat );

/*
SLA_EQGAL - J2000 alpha,delta to Galactic

ACTION: Transformation from J2000.0 FK5 equatorial coordinates to IAU 1958 galactic coordinates.
CALL:   CALL sla_EQGAL (DR, DD, DL, DB)
GIVEN:  DR,DD 	D 	J2000.0
RETURNED: DL,DB 	D 	galactic longitude and latitude (radians)

NOTE: The equatorial coordinates are J2000.0 FK5. Use the routine sla_EG50 if conversion from B1950.0 FK4 coordinates is required. 

REFERENCE:
    Blaauw et al., 1960, Mon.Not.R.astr.Soc., 121, 123. 
*/

void csla_eqgal(double dr, double dd, double *dl, double *db);

//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_EVP - Earth Position & Velocity   

ACTION:    Barycentric and heliocentric velocity and position of the Earth. 
CALL:      CALL sla_EVP (DATE, DEQX, DVB, DPB, DVH, DPH)
GIVEN:
    DATE 	D 	TDB (formerly ET) as a Modified Julian Date (JD-2400000.5)
    DEQX 	D 	Julian Epoch (e.g. 2000D0) of mean equator and equinox of the vectors returned. If DEQX <0, all vectors are referred to the mean equator and equinox (FK5) of date DATE.
RETURNED:
    DVB 	D(3) 	barycentric $[\,\dot{x},\dot{y},\dot{z}\,]$, AU s-1
    DPB 	D(3) 	barycentric $[\,x,y,z\,]$, AU
    DVH 	D(3) 	heliocentric $[\,\dot{x},\dot{y},\dot{z}\,]$, AU s-1
    DPH 	D(3) 	heliocentric $[\,x,y,z\,]$, AU

NOTES:
    1. This routine is used when accuracy is more important than CPU time, yet the extra complication of reading a pre-computed ephemeris is not justified. The maximum deviations from the JPL DE96 ephemeris are as follows:
            velocity (barycentric or heliocentric): 420 mm s-1
            position (barycentric): 6900 km
            position (heliocentric): 1600 km 
    2. The routine is an adaption of the BARVEL and BARCOR subroutines of P.Stumpff, which are described in Astr.Astrophys.Suppl.Ser. 41, 1-8 (1980). Most of the changes are merely cosmetic and do not affect the results at all. However, some adjustments have been made so as to give results that refer to the new (IAU 1976 `FK5') equinox and precession, although the differences these changes make relative to the results from Stumpff's original `FK4' version are smaller than the inherent accuracy of the algorithm. One minor shortcoming in the original routines that has not been corrected is that slightly better numerical accuracy could be achieved if the various polynomial evaluations were to be so arranged that the smallest terms were computed first. Note also that one of Stumpff's precession constants differs by 0.001" from the value given in the Explanatory Supplement. 
 */
void csla_evp(double date, double deqx, double dvb[3], double dpb[3], double dvh[3], double dph[3]);

//START REGION DEVELOP

/*
ACTION: Transformation from IAU 1958 galactic coordinates to J2000.0 FK5 equatorial coordinates.
CALL:   CALL sla_GALEQ (DL, DB, DR, DD)
GIVEN:  DL,DB 	D 	galactic longitude and latitude
RETURNED: DR,DD 	D 	J2000.0

NOTES:
    1. All arguments are in radians. 
    2. The equatorial coordinates are J2000.0 FK5. Use the routine sla_GE50 if conversion to B1950.0 FK4 coordinates is required. 
*/
void csla_galeq(double dl, double db, double *dr, double *dd);

/*
 Convert geodetic position to geocentric (double precision)

Given: P dp latitude (geodetic, radians) H dp height above reference spheroid (geodetic, metres)

Returned: R dp distance from Earth axis (AU) Z dp distance from plane of Earth equator (AU)

Notes: 1) Geocentric latitude can be obtained by evaluating ATAN2(Z,R). 2) IAU 1976 constants are used.

Reference: Green,R.M., Spherical Astronomy, CUP 1985, p98. 
*/
void csla_geoc(double p, double h, double *r, double *z);


//START REGION RELEASE

/*
SLA_GMST - UT to GMST   

ACTION: Conversion from universal time UT1 to Greenwich mean sidereal time.
CALL:   D = sla_GMST (UT1)
GIVEN:  UT1 	D 	universal time (strictly UT1) expressed as modified Julian Date (JD-2400000.5)
RETURNED: sla_GMST 	D 	Greenwich mean sidereal time (radians)

NOTES:
    1. The IAU 1982 expression (see page S15 of the 1984 Astronomical Almanac) is used, but rearranged to reduce rounding errors. This expression is always described as giving the GMST at 0h UT; in fact, it gives the difference between the GMST and the UT, which happens to equal the GMST (modulo 24 hours) at 0h UT each day. In sla_GMST, the entire UT is used directly as the argument for the canonical formula, and the fractional part of the UT is added separately; note that the factor 1.0027379... does not appear. 
    2. See also the routine sla_GMSTA, which delivers better numerical precision by accepting the UT date and time as separate arguments. 
*/
double csla_gmst(double ut1);

//START REGION DEVELOP


/*
SLA_MAP - Mean to Apparent   

ACTION:  Transform star alpha,delta from mean place to geocentric apparent. The reference frames and timescales used are post IAU 1976.
CALL:    CALL sla_MAP (RM, DM, PR, PD, PX, RV, EQ, DATE, RA, DA)
GIVEN:
    RM,DM 	D 	mean alpha,delta (radians)
    PR,PD 	D 	proper motions: alpha,delta changes per Julian year
    PX 	D 	parallax (arcsec)
    RV 	D 	radial velocity (km s-1, +ve if receding)
    EQ 	D 	epoch and equinox of star data (Julian)
    DATE 	D 	TDB for apparent place (JD-2400000.5)
RETURNED:
    RA,DA 	D 	apparent alpha,delta (radians)

NOTES:
    1.  EQ is the Julian epoch specifying both the reference frame and the epoch of the position - usually 2000. For positions where the epoch and equinox are different, use the routine sla_PM to apply proper motion corrections before using this routine. 
    2.  The distinction between the required TDB and TT is always negligible. Moreover, for all but the most critical applications UTC is adequate. 
    3.  The alpha proper motions are $\dot{\alpha}$ rather than $\dot{\alpha}\cos\delta$, and are per year rather than per century. 
    4.  This routine may be wasteful for some applications because it recomputes the Earth position/velocity and the precession/nutation matrix each time, and because it allows for parallax and proper motion. Where multiple transformations are to be carried out for one epoch, a faster method is to call the sla_MAPPA routine once and then either the sla_MAPQK routine (which includes parallax and proper motion) or sla_MAPQKZ (which assumes zero parallax and FK5 proper motion). 

REFERENCES:
    1. 1984 Astronomical Almanac, pp B39-B41. 
    2. Lederle & Schwan, 1984. Astr.Astrophys. 134, 1-6. 
 */
void csla_map(double rm, double dm, double pr, double pd, double px, double rv, double eq, double date, double *ra, double *da);

//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_NUT - Nutation Matrix   

ACTION: Form the matrix of nutation (IAU 1980 theory) for a given date. 
CALL:   CALL sla_NUT (DATE, RMATN)
GIVEN:
    DATE 	D 	TDB (formerly ET) as Modified Julian Date (JD-2400000.5)
RETURNED:
    RMATN 	D(3,3) 	nutation matrix
NOTE:
    The matrix is in the sense:
        vtrue = M.vmean
    where vtrue is the star vector relative to the true equator and equinox of date, M is the 3x3 matrix RMATN and vmean is the star vector relative to the mean equator and equinox of date. 

REFERENCES:
    1. Final report of the IAU Working Group on Nutation, chairman P.K.Seidelmann, 1980. 
    2. Kaplan, G.H., 1981. USNO circular No. 163, pA3-6. 
 */
void csla_nut(double date, double rmatn[3][3]);


//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_PA - h, delta to Parallactic Angle
ACTION: Hour angle and declination to parallactic angle (double precision).
CALL:   D = sla_PA (HA, DEC, PHI)
GIVEN:
    HA 	D 	hour angle in radians (geocentric apparent)
    DEC 	D 	declination in radians (geocentric apparent)
    PHI 	D 	latitude in radians (geodetic)
RETURNED:
    sla_PA 	D 	parallactic angle (radians, in the range +-pi)

NOTES:
    1. The parallactic angle at a point in the sky is the position angle of the vertical, i.e. the angle between the direction to the pole and to the zenith. In precise applications care must be taken only to use geocentric apparent h,delta and to consider separately the effects of atmospheric refraction and telescope mount errors. 
    2. At the pole a zero result is returned. 
*/
double csla_pa(double ha, double dec, double phi);

//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_PREBN - Precession Matrix (FK4)   

ACTION: Generate the matrix of precession between two epochs, using the old, pre IAU 1976, Bessel-Newcomb model, in Andoyer's formulation.
CALL:  CALL sla_PREBN (BEP0, BEP1, RMATP)
GIVEN:
    BEP0 	D 	beginning Besselian epoch
    BEP1 	D 	ending Besselian epoch
RETURNED:
    RMATP 	D(3,3) 	precession matrix

NOTE: The matrix is in the sense:

        v1 = M.v0

    where v1 is the star vector relative to the mean equator and equinox of epoch BEP1, M is the 3x3 matrix RMATP and v0 is the star vector relative to the mean equator and equinox of epoch BEP0. 

REFERENCE:
    Smith et al., 1989. Astr.J. 97, 269. 

 */
void csla_prebn(double bep0, double bep1, double rmatp[3][3]);

//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_PREC - Precession Matrix (FK5)   

ACTION: Form the matrix of precession between two epochs (IAU 1976, FK5). 
CALL:   CALL sla_PREC (EP0, EP1, RMATP)
GIVEN:
    EP0 	D 	beginning epoch
    EP1 	D 	ending epoch
RETURNED:
    RMATP 	D(3,3) 	precession matrix

NOTES:
    1. The epochs are TDB Julian epochs. 
    2. The matrix is in the sense:
            v1 = M.v0
        where v1 is the star vector relative to the mean equator and equinox of epoch EP1, M is the 3x3 matrix RMATP and v0 is the star vector relative to the mean equator and equinox of epoch EP0. 
    3. Though the matrix method itself is rigorous, the precession angles are expressed through canonical polynomials which are valid only for a limited time span. There are also known errors in the IAU precession rate. The absolute accuracy of the present formulation is better than 0.1" from 1960AD to 2040AD, better than 1" from 1640AD to 2360AD, and remains below 3" for the whole of the period 500BC to 3000AD. The errors exceed 10" outside the range 1200BC to 3900AD, exceed 100" outside 4200BC to 5600AD and exceed 1000" outside 6800BC to 8200AD. The SLALIB routine sla_PRECL implements a more elaborate model which is suitable for problems spanning several thousand years. 
*/
void csla_prec(double ep0, double ep1, double rmatp[3][3]);

//START REGION DEVELOP
/*
SLA_PRECES - Precession   

ACTION: Precession - either the old ``FK4'' (Bessel-Newcomb, pre IAU 1976) or new ``FK5'' (Fricke, post IAU 1976) as required.
CALL:   CALL sla_PRECES (SYSTEM, EP0, EP1, RA, DC)
GIVEN:
    SYSTEM 	C 	precession to be applied: `FK4' or `FK5'
    EP0,EP1 	D 	starting and ending epoch
    RA,DC 	D 	alpha,delta, mean equator & equinox of epoch EP0
RETURNED:
    RA,DC 	D 	alpha,delta, mean equator & equinox of epoch EP1

NOTES:
    1. Lowercase characters in SYSTEM are acceptable. 
    2. The epochs are Besselian if SYSTEM=`FK4' and Julian if `FK5'. For example, to precess coordinates in the old system from equinox 1900.0 to 1950.0 the call would be:
            CALL sla_PRECES ('FK4', 1900D0, 1950D0, RA, DC) 
    3. This routine will NOT correctly convert between the old and the new systems - for example conversion from B1950 to J2000. For these purposes see sla_FK425, sla_FK524, sla_FK45Z and sla_FK54Z. 
    4. If an invalid SYSTEM is supplied, values of -99D0,-99D0 are returned for both RA and DC. 
*/
void csla_preces(char *system, double ep0, double ep1, double *ra, double *dec);


/*
SLA_PRECL - Precession Matrix (latest)   

ACTION: Form the matrix of precession between two epochs, using the model of Simon et al. (1994), which is suitable for long periods of time.
CALL:   CALL sla_PRECL (EP0, EP1, RMATP)
GIVEN:
    EP0 	D 	beginning epoch
    EP1 	D 	ending epoch
RETURNED:
    RMATP 	D(3,3) 	precession matrix

NOTES:

    1. The epochs are TDB Julian epochs. 
    2. The matrix is in the sense:
            v1 = M.v0
        where v1 is the star vector relative to the mean equator and equinox of epoch EP1, M is the 3x3 matrix RMATP and v0 is the star vector relative to the mean equator and equinox of epoch EP0. 
    3. The absolute accuracy of the model is limited by the uncertainty in the general precession, about 0.3" per 1000 years. The remainder of the formulation provides a precision of 1 milliarcsecond over the interval from 1000AD to 3000AD, 0.1" from 1000BC to 5000AD and 1" from 4000BC to 8000AD. 

REFERENCE:
    Simon, J.L. et al., 1994. Astr.Astrophys. 282, 663. 
*/
void csla_precl(double ep0, double ep1, double rmatp[3][3]);


//START REGION DEVELOP
//START REGION RELEASE

/*
SLA_PRENUT - Precession/Nutation Matrix   

ACTION: Form the matrix of precession and nutation (IAU 1976, FK5). 
CALL:   CALL sla_PRENUT (EPOCH, DATE, RMATPN)
GIVEN:
    EPOCH 	D 	Julian Epoch for mean coordinates
    DATE 	D 	Modified Julian Date (JD-2400000.5) for true coordinates
RETURNED:
    RMATPN 	D(3,3) 	combined precession/nutation matrix

NOTES:
    1. The epoch and date are TDB. 
    2. The matrix is in the sense:
            vtrue = M.vmean
        where vtrue is the star vector relative to the true equator and equinox of epoch DATE, M is the 3x3 matrix RMATPN and vmean is the star vector relative to the mean equator and equinox of epoch EPOCH. 
 */
void csla_prenut(double ep0, double date, double rmatpn[3][3]);

//START REGION DEVELOP

/*
SLA_SUBET - Remove E-terms   

ACTION: Remove the E-terms (elliptic component of annual aberration) from a pre IAU 1976 catalogue alpha,delta to give a mean place.

CALL: CALL sla_SUBET (RC, DC, EQ, RM, DM)
GIVEN:
    RC,DC 	D 	alpha,delta with E-terms included (radians)
    EQ 	D 	Besselian epoch of mean equator and equinox
RETURNED:
    RM,DM 	D 	alpha,delta without E-terms (radians)

NOTE:
    Most star positions from pre-1984 optical catalogues (or obtained by astrometry with respect to such stars) have the E-terms built-in. This routine converts such a position to a formal mean place (allowing, for example, comparison with a pulsar timing position). 

REFERENCE:
    Explanatory Supplement to the Astronomical Ephemeris, section 2D, page 48. 
 */
void csla_subet(double rc, double dc, double eq, double *rm, double *dm);

//START REGION DEVELOP
