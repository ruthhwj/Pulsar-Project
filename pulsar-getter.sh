#!/bin/bash
# Get pulsars 
fakePulsar -v  -cone "$1 $2 $3 $4 750 12 -1 300 0" -ellipse "$5 $6 1" -cone "$7 $8 $9 ${10} 750 12 -1 300 0" -ellipse "$5 $6 1" -a "${11}" -b "${12}" -nb "2246" -np "620" -gg ${13}
pmod -templatedata weak.all37.p3fold -align -oformat ASCII -ext p3fold.ASCII ${13} 
