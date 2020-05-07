#!/bin/bash
# Get pulsars 
fakePulsar -normdisable -cone "$1 $2 $3 $4 750 $5 -1 300 0" -ellipse "$6 $7 1" -cone "$8 $9 ${10} ${11} 750 ${12} -1 300 0" -ellipse "$6 $7 1" -a "20" -b "-2.4" -nb "1123" -np "50" -gg ${13}
pmod -templatedata weak.all37.p3fold.rebinned.ASCII -align -oformat ASCII ${13} -output ${13}.ASCII
#(creates ${13}.final.ASCII)