      SUBROUTINE sla_SUBET (RC, DC, EQ, RM, DM)
*+
*     - - - - - -
*      S U B E T
*     - - - - - -
*
*  Remove the E-terms (elliptic component of annual aberration)
*  from a pre IAU 1976 catalogue RA,Dec to give a mean place
*  (double precision)
*
*  Given:
*     RC,DC     dp     RA,Dec (radians) with E-terms included
*     EQ        dp     Besselian epoch of mean equator and equinox
*
*  Returned:
*     RM,DM     dp     RA,Dec (radians) without E-terms
*
*  Called:
*     sla_ETRMS, sla_DCS2C, sla_,DVDV, sla_DCC2S, sla_DRANRM
*
*  Explanation:
*     Most star positions from pre-1984 optical catalogues (or
*     derived from astrometry using such stars) embody the
*     E-terms.  This routine converts such a position to a
*     formal mean place (allowing, for example, comparison with a
*     pulsar timing position).
*
*  Reference:
*     Explanatory Supplement to the Astronomical Ephemeris,
*     section 2D, page 48.
*
*  P.T.Wallace   Starlink   10 May 1990
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*-

      IMPLICIT NONE

      DOUBLE PRECISION RC,DC,EQ,RM,DM

      DOUBLE PRECISION sla_DRANRM,sla_DVDV
      DOUBLE PRECISION A(3),V(3),F

      INTEGER I



*  E-terms
      CALL sla_ETRMS(EQ,A)

*  Spherical to Cartesian
      CALL sla_DCS2C(RC,DC,V)

*  Include the E-terms
      F=1D0+sla_DVDV(V,A)
      DO I=1,3
         V(I)=F*V(I)-A(I)
      END DO

*  Cartesian to spherical
      CALL sla_DCC2S(V,RM,DM)

*  Bring RA into conventional range
      RM=sla_DRANRM(RM)

      END
