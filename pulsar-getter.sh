#!/bin/bash
# Get pulsars 
fakePulsar -v  -normdisable -P 1 -cone "$1 $2 $3 $4 750 17 -1 1000 0" -ellipse "$5 $6 1" -cone "$7 $8 $9 ${10} 750 29 -1 1000 0" -ellipse "$5 $6 1" -a "${11}" -b "${12}" -nb "2246" -np "50" -gg ${13}
#pulsar_arg=["./pulsar-getter.sh", "227.7", "10.5", "1","15","0.8","45","69.95","7.7","0.85", "15", "20", "-2.4", "refpulsar.gg"]
pmod -templatedata weak.all37.p3fold.ASCII -align -oformat ASCII ${13} -output ${13}.ASCII
#(creates ${13}.final.ASCII)