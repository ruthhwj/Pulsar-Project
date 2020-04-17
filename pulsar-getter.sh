#!/bin/bash
# Get pulsars 
fakePulsar -v  -normdisable -P 1 -cone "$1 $2 $3 $4 750 17 -1 1000 0" -ellipse "$5 $6 1" -cone "$7 $8 $9 ${10} 750 29 -1 1000 0" -ellipse "$5 $6 1" -a "${11}" -b "${12}" -nb "2246" -np "50" -gg ${13}

pmod -templatedata weak.all37.p3fold.ASCII -align -oformat ASCII ${13} -output ${13}.ASCII
#(creates ${13}.final.ASCII)