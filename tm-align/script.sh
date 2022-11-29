#!/usr/bin/bash

ALL_PAIRS=$1
RESULTS="$ALL_PAIRS.result"
# cat ${ALL_PAIRS} | xargs -n2 TMalign | grep -o -E "TM-score= ([+-]?[0-9]*[.]?[0-9]+)"
while read line
do
    echo $line >> $RESULTS
    TMalign $line | grep -o -E "TM-score= ([+-]?[0-9]*[.]?[0-9]+)" >> $RESULTS
    echo >> $RESULTS
done < $ALL_PAIRS
# join -j 2 -o 1.1,2.1 ${ALL_PDBS} ${ALL_PDBS} | awk '!seen[$1,$2]++ && !seen[$2,$1]++' > ${ALL_PAIRS}



