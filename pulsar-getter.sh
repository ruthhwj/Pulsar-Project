#!/bin/bash
# Get pulsars 
fakePulsar -v -normdisable -cone "$1 $2 $3 $4 750 $5 -1 300 0" -ellipse "$6 $7 1" -cone "$8 $9 ${10} ${11} 750 ${12} -1 300 0" -ellipse "$6 $7 1" -a "${13}" -b "${14}" -nb "2246" -np "50" -gg ${15}
pmod -templatedata weak.all37.p3fold.ASCII -align -oformat ASCII ${15} -output ${15}.ASCII
#(creates ${13}.final.ASCII)