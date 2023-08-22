#!/bin/csh
#Usage:
# cutslice.sh <restart_file_name> <a> <b> <c> <d0> <d1> 
# e.g. cutslice.sh rs0660 1 0 0 0 4000

#Print out usage information
if ($#argv != 6) then
    echo Usage:
    echo "  " $0 \<restart_file_name\> \<a\> \<b\> \<c\> \<d0\> \<d1\>
    echo "e.g. cutslice.sh rs0660 1 0 0 0 4000 "
    exit
endif

#main script begins
sed 's/,/ , /g' $1.pov | sed 's/</ < /g' | sed 's/>/ > /g' > $1-ext.pov

awk -v a=$2 -v b=$3 -v c=$4 -v d0=$5 -v d1=$6 -f cutslice.awk $1-ext.pov > $1-sel.pov

rm -f $1-ext.pov

