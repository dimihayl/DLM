#!/bin/bash
#open sub-shell to ensure that you do not change the default environment
(
input="./CMakeDLM.txt"
index=0

while read line ; do
    MYARRAY[$index]="$line"
    index=$(($index+1))
done < $input
THISROOT=${MYARRAY[2]}/bin/thisroot.sh
. $THISROOT
if [ "$#" -eq 0 ]; then
${MYARRAY[6]}/${MYARRAY[0]}
else
${MYARRAY[6]}/${MYARRAY[0]} $*
fi

)
