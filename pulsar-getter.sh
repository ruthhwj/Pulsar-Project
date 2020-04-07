#!/bin/bash
# Get pulsars 
fakePulsar -v  -cone "$1 $2 $3 $4 750 12 -1 300 0" -ellipse "$5 $6 1" -cone "$7 $8 $9 ${10} 750 12 -1 300 0" -ellipse "$5 $6 1" -a "${11}" -b "${12}" -nb "2246" -np "620" -gg ${13}
pmod -norm_global -dev /NULL -ext gg.D.normalised -debase ${13}
#(creates ${13}.D.normalised)
pmod -addnoise 0.0432210225 ${13}.D.normalised -output ${13}.noise
#(creates ${13}.noise)
pmod -norm_global -dev /NULL -templatedata weak.all37.p3fold.ASCII -align -oformat ASCII ${13}.noise -output ${13}.final.ASCII
#(creates ${13}.final.ASCII)